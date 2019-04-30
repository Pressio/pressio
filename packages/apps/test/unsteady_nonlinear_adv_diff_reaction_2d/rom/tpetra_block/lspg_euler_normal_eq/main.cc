
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"
#include "utils_tpetra.hpp"
#include "../../../fom/gold_states_implicit.hpp"

using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReac2dBlockTpetra;
using scalar_t		= typename app_t::scalar_type;
using uint_t		= unsigned int;
using tcomm_t		= Teuchos::MpiComm<int>;
using rcpcomm_t		= Teuchos::RCP<const tcomm_t>;
using eig_dyn_mat	= Eigen::MatrixXd;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;

constexpr double eps	= 1e-12;
constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
constexpr auto zero	= ::rompp::core::constants::zero<scalar_t>();
constexpr auto t0	= zero;


struct LSPGRunner{
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;
  using fom_state_w_t	= rompp::core::Vector<app_state_t>;
  using fom_res_w_t	= rompp::core::Vector<app_residual_t>;

  rcpcomm_t comm_;
  const int Nx_ = {};
  const int Ny_ = {};
  const int romSize_ = {};

  LSPGRunner(rcpcomm_t comm,
	     int Nx, int Ny,
	     int romSize)
    : comm_{comm},
      Nx_{Nx}, Ny_{Ny},
      romSize_{romSize}{}

  fom_state_w_t run(scalar_t dt, uint_t Nsteps)
  {
    using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;
    using mv_t		= Tpetra::Experimental::BlockMultiVector<>;
    using decoder_jac_t	= rompp::core::MultiVector<mv_t>;
    using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

    // app object
    app_t appobj(comm_, Nx_, Ny_);
    appobj.setup();
    // get the mesh map
    auto gridMap = appobj.getDataMap();
    auto ptMap = appobj.getPointMap();

    // total number of dofs (note again that this != mesh pts)
    const auto totDofs = appobj.getUnknownCount();

    // number of species
    auto const numSpecies = appobj.getNumSpecies();

    // store modes computed before from file: for now, the basis is
    // stored into a regular tpetra::MultiVector
    const auto phi0 =
      rompp::apps::test::tpetra::readBasis("basis.txt", romSize_,
					   totDofs, comm_, ptMap);

    // we now convert this MV into a Tpetra::BlockMultiVector
    mv_t phi1(*phi0.data(), *gridMap, numSpecies);
    // wrap into our MV class
    ::rompp::core::MultiVector<mv_t> phi(phi1);

    // decoder object
    decoder_t decoderObj(phi);

    // my reference state
    auto yRef = appobj.getInitialState();
    yRef.putScalar(zero);

    // define ROM state and set to zero
    lspg_state_t yROM(romSize_);
    yROM.putScalar(0.0);

    // define LSPG type
    using lspg_problem_type = rompp::rom::DefaultLSPGTypeGenerator<
      app_t, ode_case, decoder_t, lspg_state_t>;
    using lspg_generator = rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_type>;
    lspg_generator lspgProblem(appobj, yRef, decoderObj, yROM, t0);

    // linear solver
    using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t  = rompp::core::Matrix<eig_dyn_mat>;
    using solver_tag   = rompp::solvers::linear::iterative::Bicgstab;
    using linear_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using lspg_stepper_t = typename lspg_problem_type::lspg_stepper_t;
    using gnsolver_t   = rompp::solvers::iterative::GaussNewton<
      lspg_stepper_t, linear_solver_t>;
    gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(100);

    // integrate in time
    rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM,
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

  // scope guard needed (MPI init within trilinos)
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));

    // romSize: because we are doing reproducing ROM
    const int romSize = Nsteps+1;
    // object running the rom
    LSPGRunner rom(Comm, Nx, Ny, romSize);
    // run LSPG
    auto romY = rom.run(dt, Nsteps);
    auto romY_vv = romY.data()->getVectorView();
    auto romY_vvd = romY_vv.getData();

    {
      // check that solution is right
      using namespace rompp::apps::test;
      const auto trueY
	= NonLinAdvDiffReac2dImpGoldStates<ode_case>::get(Nx, Ny, dt, fint);

      if (trueY.empty()) {
	std::cout << " true solution not found, empty " << std::endl;
	checkStr = "FAILED";
      }
      for (size_t i=0; i<trueY.size(); i++){
	std::cout << std::fixed << std::setprecision(15)
		  << "i=" << i
		  << " true= " << trueY[i]
		  << " romY= " << romY_vvd[i] << std::endl;

	if (std::abs(romY_vvd[i] - trueY[i]) > eps or std::isnan(romY_vvd[i])){
	  checkStr = "FAILED";
	  break;
	}
      }
    }

  }//scope guard

  std::cout << checkStr << std::endl;
  return 0;
}
