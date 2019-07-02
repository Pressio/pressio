
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"
#include "utils_epetra.hpp"
#include "../../../fom/gold_states_implicit.hpp"

using app_t		= pressio::apps::UnsteadyNonLinAdvDiffReac2dEpetra;
using scalar_t		= typename app_t::scalar_type;
using uint_t		= unsigned int;
constexpr double eps	= 1e-12;
constexpr auto ode_case = pressio::ode::ImplicitEnum::Euler;
constexpr auto zero	= ::pressio::utils::constants::zero<scalar_t>();
constexpr auto t0	= zero;


struct LSPGRunner{
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using ode_state_t	= pressio::containers::Vector<app_state_t>;
  using ode_res_t		= pressio::containers::Vector<app_rhs_t>;
  using eig_dyn_mat	= Eigen::MatrixXd;

  const Epetra_MpiComm & comm_;
  const int Nx_ = {};
  const int Ny_ = {};
  const int romSize_ = {};

  LSPGRunner(const Epetra_MpiComm & comm,
	     int Nx, int Ny,
	     int romSize)
    : comm_{comm},
      Nx_{Nx}, Ny_{Ny},
      romSize_{romSize}{}

  ode_state_t run(scalar_t dt, uint_t Nsteps)
  {
    using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
    using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
    using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
    using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

    // app object
    app_t appobj(comm_, Nx_, Ny_);
    appobj.setup();
    const auto totDofs = appobj.getUnknownCount();

    // store modes computed before from file
    const decoder_jac_t phi =
      pressio::apps::test::epetra::readBasis("basis.txt", romSize_,
					   totDofs,
					   comm_,
					   appobj.getDataMap());

    // decoder object
    decoder_t decoderObj(phi);

    // my reference state
    auto yRef = appobj.getInitialState();
    yRef.PutScalar(zero);

    // define ROM state and set to zero
    lspg_state_t yROM(romSize_);
    yROM.putScalar(0.0);

    // define LSPG type
    using lspg_problem_type = pressio::rom::DefaultLSPGTypeGenerator<
      app_t, ode_case, decoder_t, lspg_state_t>;
    using lspg_generator = pressio::rom::LSPGUnsteadyProblemGenerator<lspg_problem_type>;
    lspg_generator lspgProblem(appobj, yRef, decoderObj, yROM, t0);

    // linear solver
    using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t  = pressio::containers::Matrix<eig_dyn_mat>;
    using solver_tag   = pressio::solvers::linear::iterative::Bicgstab;//LSCG;
    using linear_solver_t = pressio::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using lspg_stepper_t = typename lspg_problem_type::lspg_stepper_t;
    using gnsolver_t   = pressio::solvers::iterative::GaussNewton<
      lspg_stepper_t, linear_solver_t>;
    gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(100);

    // integrate in time
    pressio::ode::integrateNSteps(lspgProblem.stepperObj_, yROM,
				t0, dt, Nsteps, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);

    return yFomFinal;
  }//run
};

int main(int argc, char *argv[]){

  std::string checkStr {"PASSED"};

  // fix parameters to use for runs
  constexpr int Nx = 11, Ny = Nx*2-1;

  // for time integration
  constexpr scalar_t dt = 0.1;
  constexpr auto Nsteps = static_cast<unsigned int>(10);
  constexpr scalar_t fint = Nsteps*dt;

  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  if (Comm.NumProc() != 1){
    std::cout << "you can only run this with 1 rank" << std::endl;
    exit(EXIT_FAILURE);
  }

  // romSize: because we are doing reproducing ROM
  const int romSize = Nsteps+1;
  // object running the rom
  LSPGRunner rom(Comm, Nx, Ny, romSize);
  // run LSPG
  auto romY = rom.run(dt, Nsteps);

  {
    // check that solution is right
    using namespace pressio::apps::test;
    const auto trueY
      = NonLinAdvDiffReac2dImpGoldStates<ode_case>::get(Nx, Ny, dt, fint);

    if (trueY.empty()) {
      std::cout << " true solution not found, empty " << std::endl;
      checkStr = "FAILED";
    }
    if (trueY.size() != (size_t) romY.globalSize()) {
      std::cout << " size of true sol != size of computed one " << std::endl;
      checkStr = "FAILED";
    }
    for (size_t i=0; i<trueY.size(); i++){
      std::cout << std::fixed << std::setprecision(15)
  		<< "i=" << i
  		<< " true= " << trueY[i]
  		<< " romY= " << romY[i] << std::endl;

      if (std::abs(romY[i] - trueY[i]) > eps or std::isnan(romY[i])){
  	checkStr = "FAILED";
  	break;
      }
    }
  }

  std::cout << checkStr << std::endl;
  MPI_Finalize();
  return 0;
}
