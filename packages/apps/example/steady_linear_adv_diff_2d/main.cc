
#include "CORE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_STEADY"
#include "APPS_STEADYLINADVDIFF2D"
#include "utils_epetra.hpp"
#include <random>

using gen_t	   = std::mt19937;
using rand_distr_t = std::uniform_real_distribution<double>;
using rand_int_t   = std::uniform_int_distribution<int>;

// range of Prandtl number
constexpr std::array<double,2> Pr_range{1.0, 5.0};

// range of Reynolds number
constexpr std::array<double,2> Re_range{10., 100.0};


struct ResidualSampler{
  using vec_t = rompp::core::Vector<Epetra_Vector>;

  // max value of the index for sampling residual vector
  const int maxIndex_ = {};

  // contains 50 random indices for residual
  std::vector<int> randomIndices_;

  // vector to store the residual samples
  mutable std::vector<double> myR_;


  ResidualSampler(const int fullDof,
		  const int numSamples)
    : maxIndex_{fullDof-1},
      randomIndices_(numSamples, 0),
      myR_(numSamples, 0.0){}

  void setContainerElementsZero(){
    std::fill(myR_.begin(), myR_.end(), 0.0);
  }

  const std::vector<double> * getRPtr() const{
    return &myR_;
  };

  void generateIndices(){
    // random number generator (seeded)
    unsigned int seed = 842343;
    std::mt19937 engine(seed);
    rand_int_t distr(0., maxIndex_);

    auto gen = [&distr, &engine](){
		 return distr(engine);
	       };
    std::generate(randomIndices_.begin(), randomIndices_.end(), gen);
  }

  void observeResidualWhenSolverConverged(const vec_t & R) const{
    auto map = R.data()->Map();
    auto N = map.NumMyElements();
    assert( N == maxIndex_+1 );

    int j = 0;
    for (auto it=randomIndices_.begin(); it!=randomIndices_.end(); it++){
      myR_[j] = R[*it];
      j++;
    }
  }
};


int main(int argc, char *argv[]){
  using true_fom_t	= rompp::apps::SteadyLinAdvDiff2dEpetra;
  using fom_adapter_t	= rompp::apps::SteadyLinAdvDiff2dEpetraRomAdapter;
  using scalar_t	= typename fom_adapter_t::scalar_type;
  using native_state	= typename fom_adapter_t::state_type;

  std::string caseString;
  int nSamples;
  std::cout << "Enter <train> or <test>: ";
  std::cin >> caseString;
  std::cout << "You selected " << caseString << std::endl;
  assert( caseString == "train" or caseString == "test");

  std::cout << "Enter number of samples to run: ";
  std::cin >> nSamples;
  std::cout << "You want " << nSamples << " samples" << std::endl;
  assert( nSamples > 0 );

  // MPI
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 1);

  // the discretization to use for solver
  const int Nx = 11, Ny = Nx*2-1;
  // tot num of dof
  const int numDof = (Nx-2)*Ny;

  // number of samples to take from residual
  const int numResSamples = 50;

  // index of state to print as QoI
  constexpr int indexQoI = 33;

  // open files for MLEM python scripts
  // training_features.dat
  std::ofstream file_features;
  file_features.open( caseString + "_features.dat" );

  std::ofstream file_qoi_norm_fom;
  std::ofstream file_qoi_norm_rom;
  std::ofstream file_qoi_point_value_fom;
  std::ofstream file_qoi_point_value_rom;
  file_qoi_norm_fom.open( caseString + "_qoi_norm_fom.dat" );
  file_qoi_norm_rom.open( caseString + "_qoi_norm_rom.dat" );
  file_qoi_point_value_fom.open( caseString + "_qoi_point_value_fom.dat" );
  file_qoi_point_value_rom.open( caseString + "_qoi_point_value_rom.dat" );

  //---------------------------------------
  // generate random samples of Parameters

  // random number generator (seeded)
  // different seed if training or test
  unsigned int seed = (caseString.compare("train") == 0) ? 454438 : 279283;

  std::mt19937 engine(seed);
  rand_distr_t distr(0., 1.0);
  auto genPr = [&distr, &engine](){
		 auto c1 = Pr_range[1]-Pr_range[0];
		 auto c2 = Pr_range[0];
		 return c1 * distr(engine) + c2;
	     };

  auto genRe = [&distr, &engine](){
		 auto c1 = Re_range[1]-Re_range[0];
		 auto c2 = Re_range[0];
		 return c1 * distr(engine) + c2;
	     };

  // fill for Prandtl
  std::vector<double> PrS(nSamples);
  std::generate(PrS.begin(), PrS.end(), genPr);

  // fill for Reynolds
  std::vector<double> ReS(nSamples);
  std::generate(ReS.begin(), ReS.end(), genRe);

  if(rank==0){
    auto it1 = PrS.begin();
    auto it2 = ReS.begin();
    for( ;it2<ReS.end(); it1++, it2++)
      std::cout << std::setprecision(15)
		<< *it1 << " " << *it2
		<< "\n";
  }

  // create sampler for residual (this is passed to solver)
  using observer_t	 = ResidualSampler;
  observer_t myResidSampler(numDof, numResSamples);
  myResidSampler.generateIndices();

  // over parameter pairs, solve ROM and FOM, output and save residual
  for (int iSample=0; iSample<nSamples; iSample++)
  {
    // we run the FOM and LSPG for same values of parameters
    scalar_t Pr = PrS[iSample];
    scalar_t Re = ReS[iSample];

    // -------------------------
    // run FOM model first
    // -------------------------
    true_fom_t  appObj(Comm, Nx, Ny, Pr, Re);
    appObj.setup();
    appObj.assembleMatrix();
    appObj.fillRhs();
    appObj.solve();
    //appObj.printStateToFile("fom.txt");
    rompp::core::Vector<native_state> yFom(*appObj.getState());

    // -------------------------
    // do LSPG ROM
    // -------------------------
    //using native_state	= typename fom_adapter_t::state_type;
    using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
    using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;
    using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
    using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

    constexpr int romSize = 5;

    // app object for running rom
    fom_adapter_t  appObjROM(Comm, Nx, Ny, Pr, Re);

    // store modes from file
    const decoder_jac_t phi =
      rompp::apps::test::epetra::readBasis("basis.txt", romSize, numDof,
					   Comm, appObjROM.getDataMap());
    // decoder object
    decoder_t decoderObj(phi);

    // my reference state
    auto yRef = appObjROM.getState();
    yRef->PutScalar( static_cast<scalar_t>(0) );

    // define ROM state and set to zero
    lspg_state_t yROM(romSize);
    yROM.putScalar(0.0);

    // define LSPG type
    using lspg_problem_type = rompp::rom::DefaultLSPGSteadyTypeGenerator<
      fom_adapter_t, decoder_t, lspg_state_t>;
    using lspg_generator_t = rompp::rom::LSPGSteadyProblemGenerator<lspg_problem_type>;
    lspg_generator_t lspgProblem(appObjROM, *yRef, decoderObj, yROM);

    using lspg_stepper_t = typename lspg_problem_type::lspg_system_t;

    // linear solver
    using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t  = rompp::core::Matrix<eig_dyn_mat>;
    using solver_tag   = rompp::solvers::linear::iterative::LSCG;
    using linear_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using gnsolver_t   = rompp::solvers::iterative::GaussNewton<
      lspg_stepper_t, linear_solver_t, observer_t>;
    gnsolver_t solver(lspgProblem.systemObj_, yROM, linSolverObj, myResidSampler);
    solver.setTolerance(1e-14);
    solver.setMaxIterations(200);
    solver.solve(lspgProblem.systemObj_, yROM);

    // Write out training data needed for MLEM
    // Features: parameter values + 50 residual samples
    // retrieve the residual that was stored by the sampler
    auto storedR = myResidSampler.getRPtr();
    // write to file
    file_features << std::scientific << std::setprecision(14) << Re << " " << Pr << " " ;
    for(auto it=storedR->begin(); it !=storedR->end(); it++){
      file_features << std::scientific << std::setprecision(14) << *it << " ";
    }
    file_features << std::endl;

    // ROM quantities of interest
    auto yFomApprox = lspgProblem.yFomReconstructor_(yROM);
    const auto energyROM = rompp::core::ops::norm2(yFomApprox);
    file_qoi_norm_rom	     << std::setprecision(13) << energyROM << std::endl;
    file_qoi_point_value_rom << std::setprecision(13) << yFomApprox[indexQoI] << std::endl;

    // FOM quantities of interest
    const auto energyFOM = rompp::core::ops::norm2(yFom);
    file_qoi_norm_fom	     << std::setprecision(13) << energyFOM << std::endl;
    file_qoi_point_value_fom << std::setprecision(13) << yFom[indexQoI] << std::endl;

    // zero out elements of residual container
    myResidSampler.setContainerElementsZero();

  }//loop over params

  // close files
  file_features.close();
  file_qoi_norm_fom.close();
  file_qoi_norm_rom.close();
  file_qoi_point_value_fom.close();
  file_qoi_point_value_rom.close();

  MPI_Finalize();
  return 0;
}
