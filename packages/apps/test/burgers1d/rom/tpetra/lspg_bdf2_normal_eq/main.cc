
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_BURGERS1D"
#include "utils_tpetra.hpp"
#include "../../../fom/gold_states_implicit.hpp"

int main(int argc, char *argv[]){
  using fom_t		= rompp::apps::Burgers1dTpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;

  using decoder_jac_t	= rompp::core::MultiVector<Tpetra::MultiVector<>>;
  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
  //assert(Comm->getSize() == 2);

  // scope guard needed for tpetra
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    // app object
    constexpr int numCell = 20;
    fom_t appobj( {5.0, 0.02, 0.02}, numCell, Comm);
    appobj.setup();
    auto t0 = static_cast<scalar_t>(0);
    scalar_t dt = 0.01;

    // read from file the jacobian of the decoder
    constexpr int romSize = 11;
    // store modes computed before from file
    decoder_jac_t phi =
      rompp::apps::test::tpetra::readBasis("basis.txt", romSize, numCell,
					   Comm, appobj.getDataMap());
    // create decoder obj
    decoder_t decoderObj(phi);

    // for this problem, my reference state = initial state
    auto & yRef = appobj.getInitialState();

    // define ROM state and initialize to zero (this has to be done)
    lspg_state_t yROM(romSize);
    yROM.putScalar(0.0);

    // define LSPG type
    constexpr auto ode_case = rompp::ode::ImplicitEnum::BDF2;
    using lspg_problem_types = rompp::rom::DefaultLSPGTypeGenerator<
      fom_t, ode_case, decoder_t, lspg_state_t>;
    rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_types> lspgProblem
      (appobj, yRef, decoderObj, yROM, t0);

    using rom_stepper_t = typename lspg_problem_types::rom_stepper_t;

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
      rom_stepper_t, linear_solver_t>;
    gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(200);

    // integrate in time
    rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);
    auto yFF_v = yFomFinal.data()->getData();

    // this is a reproducing ROM test, so the final reconstructed state
    // has to match the FOM solution obtained with bdf2, same time-step, for 10 steps
    int shift = (rank==0) ? 0 : 10;
    const int myn = yFomFinal.getDataMap().getNodeNumElements();
    const auto trueY = rompp::apps::test::Burgers1dImpGoldStates<ode_case>::get(numCell, dt, 0.10);
    for (auto i=0; i<myn; i++){
      assert(std::abs(yFF_v[i] - trueY[i+shift]) < 1e-10 );
    }

  }//tpetra scope

  MPI_Finalize();
  return 0;
}
