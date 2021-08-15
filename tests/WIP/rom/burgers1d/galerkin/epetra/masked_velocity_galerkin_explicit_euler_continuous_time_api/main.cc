
#include "pressio_rom_galerkin.hpp"
#include "pressio_apps.hpp"
#include "utils_epetra.hpp"

// note that the gold solution is one the full mesh
const std::vector<double> gold = {
  1.2397709171427,
  1.0046222839216,
  1.0000008784156,
  0.99999103326939,
  1.0032171867377,
  1.0018031000762,
  0.99834220537621,
  1.0043668267306,
  0.99913640885698,
  0.99761548539345,
  1.0026066654337,
  1.0081714325534,
  1.0053729703494,
  1.0029387880013,
  1.0039676445756,
  1.0076713526184,
  1.0040394420228,
  1.0113275860806,
  1.0108300744443,
  1.0131357962777};

using fom_t		= pressio::apps::Burgers1dEpetra;
using scalar_t		= typename fom_t::scalar_type;
using native_state_t	= typename fom_t::state_type;
using native_velo_t	= typename fom_t::velocity_type;
using fom_state_t	= pressio::containers::Vector<native_state_t>;
using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
using decoder_t		= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
using rom_state_t	= pressio::containers::Vector<Eigen::VectorXd>;

struct MyCollocator
{
  const Epetra_Import & importer_;
  MyCollocator(const Epetra_Import & importer) : importer_(importer){}

  Epetra_MultiVector sampleRows(const Epetra_MultiVector & operand)
  {
    Epetra_MultiVector result(importer_.TargetMap(), operand.NumVectors());
    result.Import(operand, importer_, Epetra_CombineMode::Insert);
    return result;
  }
};

struct Masker
{
  const Epetra_Import & importer_;
  Masker(const Epetra_Import & importer) : importer_(importer){}

  native_velo_t createApplyMaskResult(const native_velo_t & src) const
  {
    return native_velo_t(importer_.TargetMap());
  }

  void applyMask(const native_velo_t & src, double /*t*/, native_velo_t & dest) const
  {
    dest.Import(src, importer_, Epetra_CombineMode::Insert);
  }
};

Epetra_Map createCollocationMap(const Epetra_MpiComm & comm, int rank)
{
  /*
    the full grid is distributed over 3 ranks such as:
    - rank0 contains gIDs 0...6
    - rank1 contains gIDs 7...13
    - rank2 contains gIDs 14...19

    only pick subset of grid points as follows
    - rank0 contains gIDs 0,1,4,5
    - rank1 contains gIDs 7,10,11,12
    - rank2 contains gIDs 15,17,18,19
  */

  std::array<int,4> myGIDs;
  if (rank==0){
    myGIDs = {{0,1,4,5}};
  }
  else if (rank==1){
    myGIDs = {{7,10,11,12}};
  }
  else if (rank==2){
    myGIDs = {{15,17,18,19}};
  }
  return Epetra_Map(12, myGIDs.size(), myGIDs.data(), 0, comm);
}

int main(int argc, char *argv[])
{
  std::string checkStr {"PASSED"};

  // MPI init
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  if(Comm.NumProc() != 3) return 0;

  // app object
  int numCell = 20;
  fom_t appobj({5.0, 0.02, 0.02}, numCell, &Comm);
  scalar_t dt = 0.01;

  // read the basis
  constexpr int romSize = 11;
  decoder_jac_t phi = pressio::rom::test::epetra::readBasis
    ("basis.txt", romSize, numCell, Comm, appobj.getDataMap());
  decoder_t decoderObj(phi);

  // this is my reference state
  auto & y0n = appobj.getInitialState();

  // ROM state
  rom_state_t yROM(romSize);
  pressio::ops::set_zero(yROM);

  // create a map for doing the collocation
  const auto collMap = createCollocationMap(Comm, rank);
  const Epetra_Import importer(collMap, appobj.getDataMap());

  // since this is a masked Galerkin problem, we need the collocator to
  // create a projector suitable for the masked operator
  MyCollocator mapper(importer);
  auto projector = pressio::rom::galerkin::createCollocationProjector(decoderObj,
								      std::move(mapper));

  // create masker
  Masker myMask(importer);

  using ode_tag = pressio::ode::ForwardEuler;
  auto galerkinProb = pressio::rom::galerkin::createMaskedVelocityProblem
    <ode_tag>(appobj, decoderObj, yROM, y0n, myMask, projector);

  pressio::rom::galerkin::solveNSteps(galerkinProb, yROM, 0.0, dt, 10);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = galerkinProb.fomStateReconstructorCRef()(yROM);
  yFomFinal.data()->Print(std::cout << std::setprecision(14));

  // check against gold solution
  int shift = 0;
  if (rank==1)  shift = 7;
  if (rank==2)  shift = 14;
  int myn = yFomFinal.data()->Map().NumMyElements();
  for (auto i=0; i<myn; i++){
    if(std::abs(yFomFinal(i) - gold[i+shift]) > 1e-12 ){
      checkStr = "FAILED";
      break;
    }
  }

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;
}
