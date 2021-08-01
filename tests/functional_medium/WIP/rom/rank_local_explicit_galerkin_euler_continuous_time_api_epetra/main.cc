
#include "pressio_rom_galerkin.hpp"
#include "pressio_apps.hpp"

namespace{

class FakeApp
{
public:
  using nativeVec = Epetra_Vector;
  template<typename T> using rcp = std::shared_ptr<T>;

  using scalar_type = double;
  using state_type  = Epetra_Vector;
  using velocity_type = state_type;

public:
  FakeApp() = delete;
  FakeApp(int Ncell, Epetra_MpiComm * comm, std::string & sentinel)
    : Ncell_(Ncell), comm_(comm), sentinel_(sentinel)
  {
    dataMap_ = std::make_shared<Epetra_Map>(Ncell_, 0, *comm_);
    NumMyElem_ = dataMap_->NumMyElements();
    myRank_ =  comm_->MyPID();
    totRanks_ =  comm_->NumProc();
  }

  ~FakeApp() = default;

public:
  Epetra_Map const & getDataMap(){
    return *dataMap_;
  };

  velocity_type createVelocity() const{
    Epetra_Vector R(*dataMap_);
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

    for (int i=0; i<NumMyElem_; ++i){
      if ( std::abs(u[i] - goldU(i)) > 1e-13) sentinel_ = "FAILED";
    }
    //u.Print(std::cout);

    f.PutScalar(myRank_);
  }

private:
  int Ncell_;
  Epetra_MpiComm * comm_;
  std::string & sentinel_;
  rcp<Epetra_Map> dataMap_;
  int myRank_;
  int totRanks_;
  int NumMyElem_;
  std::vector<int> myGel_;
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

  std::string checkStr {"PASSED"};

  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  if(Comm.NumProc() != 4) return 0;

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  const int numCellGlobal = 40;
  fom_t appobj(numCellGlobal, &Comm, checkStr);

  native_state_t y0n(appobj.getDataMap());
  y0n.PutScalar(1.0);
  const auto numCellLocal = y0n.MyLength();

  constexpr int romSize = 4;
  decoder_jac_t phi(numCellLocal, romSize);
  // fill with a MPI rank-dependent value, for testing purposes
  pressio::ops::fill(phi, (double) rank+1.);
  decoder_t decoderObj(phi);

  rom_state_t yROM(romSize);
  pressio::ops::fill(yROM, 1.0);

  using ode_tag = pressio::ode::explicitmethods::Euler;
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
  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;
}


// #include "utils_epetra.hpp"
// #include "utils_eigen.hpp"

// // const std::vector<double> bdf1Sol
// // {5.0081542681376, 5.016629490569, 5.025433912557, 5.0345792953115,
// //  5.0440827179355, 5.0539568295087, 5.0642107801363, 5.074857742734,
// //  5.0859146000515, 5.0974001265619, 5.1093302710968, 5.1217197481536,
// //  5.1345846667293, 5.1479436063682, 5.1618137609004, 5.1762071980595,
// //  5.1911395190849, 5.2066322357211, 5.222706587389, 5.2393822195142,
// //  5.2566784890019, 5.274617970535, 5.2932246323729, 5.3125186218141,
// //  5.3325236627322, 5.3532729201416, 5.3747971779128, 5.3971189932731,
// //  5.4202577535351, 5.4442348269811, 5.469078757402, 5.4948202159562,
// //  5.5214859714822, 5.5491009348394, 5.5776911098501, 5.6072849195866,
// //  5.6379131952825, 5.6696069037791, 5.7023980878343, 5.7363239274031,
// //  5.7714263431002, 5.807744410524, 5.8453128737429, 5.884168702448,
// //  5.9243510856362, 5.9658923478856, 6.0088164545724, 6.0531503069487,
// //  6.0989210765093, 6.1461565470309};

// template <typename state_t>
// struct observer
// {
//   using matrix_t = Eigen::MatrixXd;

//   int snapFreq_ = {};
//   size_t local_state_size_ {};
//   matrix_t A_;
//   size_t count_ = 0;
//   Eigen::VectorXd localInitState_;

//   observer(int N, int state_size, int snapFreq)
//     : snapFreq_(snapFreq),
//       local_state_size_(state_size),
//       A_(state_size, (N/snapFreq)+1),
//       localInitState_(state_size)
//   {
//     localInitState_.setConstant(1.);
//   }

//   void operator()(size_t step,
//   		  double t,
//   		  const state_t & y)
//   {
//     if (step % snapFreq_ == 0){
//       for (std::size_t i=0; i<local_state_size_; ++i){
// 	A_(i, count_) = y(i) - localInitState_(i);
//       }
//       count_++;
//     }
//   }

//   size_t getCount() const{ return count_;}
//   void printAll() const{ std::cout << A_ << std::endl; }

//   void printToFile(int rank) const
//   {
//     std::ofstream file;
//     file.open("snap_"+std::to_string(rank)+".txt");
//     for (int i=0; i<A_.rows(); i++){
//       for (int j=0; j<A_.cols(); j++){
// 	file << std::setprecision(15) << A_(i,j) << " ";
//       }
//       file << std::endl;
//     }
//     file.close();
//   }
// };

// template<class fom_t, typename sc_t>
// void runFom(int rank, const fom_t & appobj, int snapFreq, sc_t dt, sc_t finalTime)
// {
//   using scalar_t	= typename fom_t::scalar_type;
//   using native_state_t  = typename fom_t::state_type;
//   using fom_state_t     = pressio::containers::Vector<native_state_t>;

//   const auto & y0n = appobj.getInitialState();
//   fom_state_t y(y0n);
//   auto stepperObj = pressio::ode::createForwardEulerStepper(y, appobj);

//   auto Nsteps = static_cast<::pressio::ode::step_type>(finalTime/dt);
//   observer<fom_state_t> observer(Nsteps, y.extentLocal(0), snapFreq);
//   pressio::ode::advanceNSteps(stepperObj, y, 0.0, dt, Nsteps, observer);
//   observer.printToFile(rank);
// }
