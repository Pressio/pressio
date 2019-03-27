
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_BURGERS1D"
#include "utils_epetra.hpp"
#include "../../../fom/fom_gold_states.hpp"

int main(int argc, char *argv[]){
  using fom_t		= rompp::apps::Burgers1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;

  using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);

  //-------------------------------
  // app object
  constexpr int numCell = 20;
  //std::vector<double> mu({5.0, 0.02, 0.02});
  fom_t appobj( {5.0, 0.02, 0.02}, numCell, &Comm);
  appobj.setup();
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;

  // read from file the jacobian of the decoder
  constexpr int romSize = 11;
  // store modes computed before from file
  decoder_jac_t phi =
    rompp::apps::test::epetra::readBasis("basis.txt", romSize, numCell,
					 Comm, appobj.getDataMap());
  const int numBasis = phi.globalNumVectors();
  assert( numBasis == romSize );

  // this is my reference state
  auto & y0n = appobj.getInitialState();
  decoder_t decoderObj(phi);

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  // define LSPG type
  using lspg_problem_types = rompp::rom::DefaultLSPGTypeGenerator<
    fom_t, rompp::ode::ImplicitEnum::BDF2, decoder_t, lspg_state_t>;
  rompp::rom::LSPGStepperObjectGenerator<lspg_problem_types> lspgGener(
      appobj, y0n, decoderObj, yROM, t0);

  using rom_residual_t = typename lspg_problem_types::lspg_residual_t;
  using rom_jac_t      = typename lspg_problem_types::lspg_matrix_t;

  // GaussNewton solver
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  // rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t	 = rompp::core::Matrix<eig_dyn_mat>;
  using solver_tag	 = rompp::solvers::linear::LSCG;
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using gnsolver_t	 = rompp::solvers::iterative::GaussNewton<
    scalar_t, solver_tag, rompp::solvers::EigenIterative,
    converged_when_t, void, hessian_t, lspg_state_t, rom_residual_t, rom_jac_t>;
  gnsolver_t solver(lspgGener.stepperObj_, yROM);
  solver.setTolerance(1e-13);
  solver.setMaxIterations(200);

  // integrate in time
  rompp::ode::integrateNSteps(lspgGener.stepperObj_, yROM, 0.0, dt, 10, solver);

  // compute the fom corresponding to our rom final state
  using fom_state_w_t = typename lspg_problem_types::fom_state_w_t;
  fom_state_w_t yRf(y0n);
  decoderObj.applyMapping(yROM, yRf);
  yRf += lspgGener.y0Fom_;
  yRf.data()->Print(std::cout << std::setprecision(14));

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with bdf2, same time-step, for 10 steps
  {
    int shift = (rank==0) ? 0 : 10;
    const int myn = yRf.getDataMap().NumMyElements();
    const auto trueY = rompp::apps::test::Burg1DtrueImpBDF2N20t010;
    for (auto i=0; i<myn; i++)
      assert(std::abs(yRf[i] - trueY[i+shift]) < 1e-10 );
  }

  MPI_Finalize();
  return 0;
}
