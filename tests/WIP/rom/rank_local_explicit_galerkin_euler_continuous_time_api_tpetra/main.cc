
#include "pressio_rom_galerkin.hpp"
#include "pressio_apps.hpp"
#include "utils_tpetra.hpp"

namespace{

class FakeApp
{
public:
  using map_t		= Tpetra::Map<>;
  using nativeVec	= Tpetra::Vector<>;
  using go_t		= typename map_t::global_ordinal_type;
  using lo_t		= typename map_t::local_ordinal_type;
  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  using rcpmap_t	= Teuchos::RCP<const map_t>;

  template<typename T> using stdrcp = std::shared_ptr<T>;
  using crs_graph_type = Tpetra::CrsGraph<>;

public:
  using scalar_type	= double;
  using state_type	= nativeVec;
  using velocity_type	= state_type;
  using jacobian_type	= Tpetra::CrsMatrix<>;

public:
  FakeApp() = delete;
  FakeApp(int Ncell, rcpcomm_t comm, std::string & sentinel)
    : Ncell_(Ncell), comm_(comm), sentinel_(sentinel)
  {
    myRank_ =  comm_->getRank();
    totRanks_ =  comm_->getSize();
    dataMap_ = Teuchos::rcp(new map_t(Ncell_, 0, comm_));
    NumMyElem_ = dataMap_->getNodeNumElements();
  }

  ~FakeApp() = default;

public:
  rcpmap_t getDataMap(){
    return dataMap_;
  };

  velocity_type createVelocity() const{
    velocity_type R(dataMap_);
    return R;
  }

  void velocity(const state_type & u,
		const scalar_type /* t */,
		velocity_type & f) const
  {
    Eigen::VectorXd goldU(NumMyElem_);
    if (myRank_==0){
      goldU.setConstant(5.);
    }
    if (myRank_==1){
      goldU.setConstant(9.);
    }
    if (myRank_==2){
      goldU.setConstant(13.);
    }
    if (myRank_==3){
      goldU.setConstant(17.);
    }

    auto uvH = u.getLocalViewHost();
    for (int i=0; i<NumMyElem_; ++i){
      if ( std::abs(uvH(i,0) - goldU(i)) > 1e-13) sentinel_ = "FAILED";
    }

    f.putScalar(myRank_);
  }

private:
  int Ncell_;
  rcpcomm_t comm_{};
  std::string & sentinel_;
  rcpmap_t dataMap_{};
  int myRank_;
  int totRanks_;
  int NumMyElem_;
  std::vector<go_t> myGel_;
};
}//end anonym namespace

int main(int argc, char *argv[])
{

  using fom_t		= FakeApp;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t     = pressio::containers::Vector<native_state_t>;
  using eig_dyn_vec     = Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t     = pressio::containers::Vector<eig_dyn_vec>;

  // note that the basis type is Eigen matrix wrapper
  // because the pod modes are local to each MPI rank
  // it is sharedmem and it is compatbile with romstate
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;

  std::string checkStr {"PASSED"};

  // scope guard needed for tpetra
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));

    pressio::log::initialize(pressio::logto::terminal);
    pressio::log::setVerbosity({pressio::log::level::debug});

    const int numCellGlobal = 40;
    fom_t appobj(numCellGlobal, Comm, checkStr);

    native_state_t y0n(appobj.getDataMap());
    y0n.putScalar(1.0);
    const auto numCellLocal = y0n.getLocalLength();

    constexpr int romSize = 4;
    decoder_jac_t phi(numCellLocal, romSize);
    // fill with a MPI rank-dependent value, for testing purposes
    pressio::ops::fill(phi, (double) rank+1.);
    decoder_t decoderObj(phi);

    rom_state_t yROM(romSize);
    pressio::ops::fill(yROM, 1.0);

    using ode_tag = pressio::ode::ForwardEuler;
    auto galerkinProb =
      pressio::rom::galerkin::createDefaultProblem<ode_tag>(appobj,
							  decoderObj,
							  yROM, y0n);

    pressio::rom::galerkin::solveNSteps(galerkinProb, yROM, 0.0, 2., 1);

    std::cout << rank << " " << *yROM.data() << std::endl;
    Eigen::VectorXd goldRom(romSize);
    if (rank==0){
      goldRom.setConstant(1.);
    }
    if (rank==1){
      goldRom.setConstant(41.);
    }
    if (rank==2){
      goldRom.setConstant(121.);
    }
    if (rank==3){
      goldRom.setConstant(241.);
    }

    if (!yROM.data()->isApprox(goldRom)){
      checkStr = "FAILED";
    }
    pressio::log::finalize();
  }

  std::cout << checkStr <<  std::endl;
  return 0;
}
