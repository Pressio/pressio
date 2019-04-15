
#include "CORE_ALL"
#include "ODE_ALL"
#include "ROM_LSPG"
#include "SOLVERS_NONLINEAR"
#include "APPS_UNSTEADYLINADVDIFF1D"
#include "utils_epetra.hpp"
#include "../fom_gold_states_implicit_euler.hpp"

constexpr double eps = 1e-7;

template <typename T>
void checkSol(int rank, const T &y,
	      const std::vector<double> & trueS){
  if (rank == 0){
    assert(std::abs(y[0] - trueS[0]) < eps);
    assert(std::abs(y[1] - trueS[1]) < eps);
    assert(std::abs(y[2] - trueS[2]) < eps);
    assert(std::abs(y[3] - trueS[3]) < eps);
    assert(std::abs(y[4] - trueS[4]) < eps);
    assert(std::abs(y[5] - trueS[5]) < eps);
    assert(std::abs(y[6] - trueS[6]) < eps);
    assert(std::abs(y[7] - trueS[7]) < eps);
    assert(std::abs(y[8] - trueS[8]) < eps);
    assert(std::abs(y[9] - trueS[9]) < eps);
    assert(std::abs(y[10] - trueS[10]) < eps);
    assert(std::abs(y[11] - trueS[11]) < eps);
    assert(std::abs(y[12] - trueS[12]) < eps);
    assert(std::abs(y[13] - trueS[13]) < eps);
    assert(std::abs(y[14] - trueS[14]) < eps);
    assert(std::abs(y[15] - trueS[15]) < eps);
    assert(std::abs(y[16] - trueS[16]) < eps);
    assert(std::abs(y[17] - trueS[17]) < eps);
    assert(std::abs(y[18] - trueS[18]) < eps);
    assert(std::abs(y[19] - trueS[19]) < eps);
  }
  if (rank == 1){
    assert(std::abs(y[0] - trueS[20]) < eps);
    assert(std::abs(y[1] - trueS[21]) < eps);
    assert(std::abs(y[2] - trueS[22]) < eps);
    assert(std::abs(y[3] - trueS[23]) < eps);
    assert(std::abs(y[4] - trueS[24]) < eps);
    assert(std::abs(y[5] - trueS[25]) < eps);
    assert(std::abs(y[6] - trueS[26]) < eps);
    assert(std::abs(y[7] - trueS[27]) < eps);
    assert(std::abs(y[8] - trueS[28]) < eps);
    assert(std::abs(y[9] - trueS[29]) < eps);
    assert(std::abs(y[10] - trueS[30]) < eps);
    assert(std::abs(y[11] - trueS[31]) < eps);
    assert(std::abs(y[12] - trueS[32]) < eps);
    assert(std::abs(y[13] - trueS[33]) < eps);
    assert(std::abs(y[14] - trueS[34]) < eps);
    assert(std::abs(y[15] - trueS[35]) < eps);
    assert(std::abs(y[16] - trueS[36]) < eps);
    assert(std::abs(y[17] - trueS[37]) < eps);
    assert(std::abs(y[18] - trueS[38]) < eps);
  }

}

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
  //----------------------------------------------------------------------
  // MPI initialization
  //----------------------------------------------------------------------
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);

  //----------------------------------------------------------------------
  // Designate the paramerters and problem
  //----------------------------------------------------------------------
  scalar_t dt = 0.1;
  scalar_t fint = 5.0;
  // Actual parameter inputs that were used to formulate the basis
  std::vector<scalar_t> mu{-0.857241631161166, 0.104833925269630,
      -0.713183149274631};
  std::vector<scalar_t> domain{0, 2.0, 0.05};
  std::vector<scalar_t> bc1D{0, 0.25};

  //----------------------------------------------------------------------
  // App for UnsteadyLinAdvDiff1dEpetra Object
  //----------------------------------------------------------------------
  auto Nsteps = static_cast<unsigned int>(fint/dt);
  fom_t appObj(Comm, mu, domain, bc1D);
  appObj.unsteadySetup();

  //----------------------------------------------------------------------
  // ROM parameters and Basis functions
  //----------------------------------------------------------------------
  constexpr int romSize = 15;
  constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
  auto t0 = static_cast<scalar_t>(0);
  const int numDof = appObj.getNumGlobalNodes();
  std::cout << numDof << "x" << romSize <<std::endl;

  //Retrieve basis found from an external source (MATLAB)
  decoder_jac_t phi =
    rompp::apps::test::epetra::readBasis("basis.txt", romSize, numDof,
					 Comm, appObj.getDataMap());
  phi.data()->Print(std::cout);
  decoder_t decoderObj(phi);

  //Retrieve the initial states defined specifically for this problem
  auto yRef =  *appObj.getInitialState();
  lspg_state_t yROM(romSize);

  //initialize to zero since using y-y0 notation
  yROM.putScalar(0.0);

  //----------------------------------------------------------------------
  // Define LSPG Type
  //----------------------------------------------------------------------
  using lspg_problem_types = rompp::rom::DefaultLSPGTypeGenerator<
    fom_t, ode_case, decoder_t, lspg_state_t>;
  rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_types>
    lspgProblem(appObj, yRef, decoderObj, yROM, t0);

  //----------------------------------------------------------------------
  // LSPG Solver
  //----------------------------------------------------------------------
  //Set solver parameters
  using rom_residual_t = typename lspg_problem_types::lspg_residual_t;
  using rom_jac_t = typename lspg_problem_types::lspg_matrix_t;
  using eig_dyn_mat = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t = rompp::core::Matrix<eig_dyn_mat>;
  using solver_tag = rompp::solvers::linear::iterative::LSCG;
  using converged_when_t = rompp::solvers::iterative::default_convergence;

  //gauss-newton solver normal equations
  using gnsolver_t       = rompp::solvers::iterative::GaussNewton<
    scalar_t, solver_tag, rompp::solvers::EigenIterative,
    converged_when_t, void, hessian_t, lspg_state_t, rom_residual_t,
    rom_jac_t>;
  gnsolver_t solver(lspgProblem.stepperObj_, yROM);
  solver.setTolerance(1e-13);
  solver.setMaxIterations(200);

  // integrate in time
  rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM, 0.0, dt,
			      Nsteps, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);
  yFomFinal.data()->Print(std::cout << std::setprecision(14));

  //----------------------------------------------------------------------
  // Compare the gold states from BDF1 to the reconstructed ROM
  //----------------------------------------------------------------------
  {
    using namespace rompp::apps::test;
    checkSol(rank, yFomFinal,
	     UnsteadyLinAdvDiff1dImplicitGoldStates<ode_case>
	     ::get(mu, domain, bc1D, dt, fint));
  }
  //----------------------------------------------------------------------
  // Finalize and Return
  //----------------------------------------------------------------------
  MPI_Finalize();
  return 0;
}
