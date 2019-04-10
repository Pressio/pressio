
#include "CORE_ALL"
#include "ODE_ALL"
#include "ROM_LSPG"
#include "SOLVERS_NONLINEAR"
#include "APPS_UNSTEADYLINADVDIFF1D"
#include "utils_epetra.hpp"
//Need to add gold solution or something here

int main(int argc, char *argv[]){
  using fom_t          = rompp::apps::UnsteadyLinAdvDiff1dEpetra;
  using scalar_t       =typename fom_t::scalar_type;
  using  eig_dyn_vec   = Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t   = rompp::core::Vector<eig_dyn_vec>;
  using decoder_jac_t  = rompp::core::MultiVector<Epetra_MultiVector>;
  using decoder_t      = rompp::rom::LinearDecoder<decoder_jac_t>;
  using native_state    = typename fom_t::state_type;
  using app_state_t     = typename fom_t::state_type;
  using app_residual_t = typename fom_t::residual_type;
  using app_jacobian_t = typename fom_t::jacobian_type;
  //---------------------------------------------------------------------------
  // MPI initialization
  //---------------------------------------------------------------------------
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);
  //---------------------------------------------------------------------------
  // Designate the paramerters and problem
  //---------------------------------------------------------------------------
  scalar_t dt = 0.1;
  scalar_t fint = 5.0;
  std::vector<scalar_t> mu{-0.857241631161166, 0.104833925269630,
      -0.713183149274631};
  std::vector<scalar_t> domain{0, 2.0, 0.05};
  std::vector<scalar_t> bc1D{0, 0.25};
  //---------------------------------------------------------------------------
  // App for UnsteadyLinAdvDiff1dEpetra Object
  //---------------------------------------------------------------------------
  fom_t appObj(Comm, mu, domain, bc1D);
  appObj.setup();
  appObj.unsteadySetup();
  // //---------------------------------------------------------------------------
  // // FOM to compare ROM with
  // //---------------------------------------------------------------------------


  // const native_state y0n(*appObj.getInitialState());
  // //  const native_state y0n(*appObj.getInitialState());
  // auto r0n = appObj.residual(y0n, static_cast<scalar_t>(0));
  // auto j0n = appObj.jacobian();

  // using ode_state_t = rompp::core::Vector<app_state_t>;
  // using ode_res_t  = rompp::core::Vector<app_residual_t>;
  // using ode_jac_t = rompp::core::Matrix<app_jacobian_t>;
  // ode_state_t yFOM(y0n);
  // ode_res_t rFOM(r0n);
  // ode_jac_t jFOM(j0n);

  // using aux_stepper_t = rompp::ode::ImplicitStepper<
  //   rompp::ode::ImplicitEnum::Euler,
  //   ode_state_t, ode_res_t, ode_jac_t, fom_t>;
  // aux_stepper_t stepperAux(yFOM, appObj);

  constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;

  // using stepper_t = rompp::ode::ImplicitStepper<
  //   ode_case, ode_state_t, ode_res_t, ode_jac_t, fom_t, aux_stepper_t, void, void>;
  // stepper_t stepperObj(yFOM, appObj, stepperAux);

  // //Integrate in time
   auto Nsteps = static_cast<unsigned int>(fint/dt);
  // rompp::ode::integrateNSteps(stepperObj, yFOM, 0.0, dt, Nsteps);
  //---------------------------------------------------------------------------
  // ROM parameters and Basis functions
  //---------------------------------------------------------------------------
  constexpr int romSize = 15;
  auto t0 = static_cast<scalar_t>(0);
  const int numDof = appObj.getNumGlobalNodes();
  std::cout << numDof << "x" << romSize <<std::endl;
  decoder_jac_t phi =
    rompp::apps::test::epetra::readBasis("basis.txt", romSize, numDof, Comm,
  					 appObj.getDataMap());
  phi.data()->Print(std::cout);
  decoder_t decoderObj(phi);

  auto yRef =  *appObj.getInitialState();

  lspg_state_t yROM(romSize);
  //initialize to zero since using y-y0 notation
  yROM.putScalar(0.0);

  //---------------------------------------------------------------------------
  // Define LSPG Type
  //---------------------------------------------------------------------------
  using lspg_problem_types = rompp::rom::DefaultLSPGTypeGenerator<
    fom_t, ode_case, decoder_t, lspg_state_t>;
  rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_types>
    lspgProblem(appObj, yRef, decoderObj, yROM, t0);
  //---------------------------------------------------------------------------
  //
  //---------------------------------------------------------------------------
  using rom_residual_t = typename lspg_problem_types::lspg_residual_t;
  using rom_jac_t = typename lspg_problem_types::lspg_matrix_t;

  using eig_dyn_mat = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t = rompp::core::Matrix<eig_dyn_mat>;
  using solver_tag = rompp::solvers::linear::iterative::LSCG;
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using gnsolver_t       = rompp::solvers::iterative::GaussNewton<
    scalar_t, solver_tag, rompp::solvers::EigenIterative, converged_when_t,
    void, hessian_t, lspg_state_t, rom_residual_t, rom_jac_t>;
  gnsolver_t solver(lspgProblem.stepperObj_, yROM);
  solver.setTolerance(1e-13);
  solver.setMaxIterations(200);

  // integrate in time
  rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM, 0.0, dt, Nsteps,
  			      solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);
  yFomFinal.data()->Print(std::cout << std::setprecision(14));


  //---------------------------------------------------------------------------
  // Compare
  //---------------------------------------------------------------------------
  // auto errorVec(yFOM);
  // errorVec = yFOM-yFomFinal;
  // const auto norm2err = rompp::core::ops::norm2(errorVec);
  // assert(norm2err < 1e-5);
  // std::cout << std::setprecision(15) << norm2err << std::endl;
  MPI_Finalize();
  return 0;
}
