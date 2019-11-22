
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "QR_BASIC"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_UNSTEADY"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_epetra.hpp"

int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  std::string checkStr {"PASSED"};

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  if(Comm.NumProc() != 1) return 0;

  //-------------------------------
  // app object
  constexpr int numCell = 20;
  fom_t appobj( {5.0, 0.02, 0.02}, numCell, &Comm);
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;

  // read from file the jacobian of the decoder
  constexpr int romSize = 11;
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::rom::test::epetra::readBasis("basis.txt", romSize, numCell,
					 Comm, appobj.getDataMap());
  const int numBasis = phi.globalNumVectors();
  if( numBasis != romSize ) return 0;

  // create decoder obj
  decoder_t decoderObj(phi);

  // for this problem, my reference state = initial state
  auto & yRef = appobj.getInitialState();

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  // define LSPG type
  constexpr auto ode_case  = pressio::ode::ImplicitEnum::BDF2;
  using lspg_problem = pressio::rom::LSPGUnsteadyProblem<
    pressio::rom::DefaultLSPGUnsteady, ode_case, fom_t, lspg_state_t, decoder_t>;
  lspg_problem lspgProblem(appobj, yRef, decoderObj, yROM, t0);

  using lspg_stepper_t = typename lspg_problem::lspg_stepper_t;
  using rom_jac_t      = typename lspg_problem::lspg_matrix_t;

  // GaussNewton solver
  using qr_algo = pressio::qr::TSQR;
  using qr_type = pressio::qr::QRSolver<rom_jac_t, qr_algo>;
  using converged_when_t = pressio::solvers::iterative::default_convergence;
  using gnsolver_t  = pressio::solvers::iterative::GaussNewtonQR<
         lspg_stepper_t, qr_type, converged_when_t>;
  gnsolver_t solver(lspgProblem.getStepperRef(), yROM);
  solver.setTolerance(1e-13);
    // I know this should converge in few iters at every step
  solver.setMaxIterations(5);

  // integrate in time
  pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM, 0.0, dt, 10, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM);
  yFomFinal.data()->Print(std::cout << std::setprecision(14));

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with bdf2, same time-step, for 10 steps
  {
    int shift = (rank==0) ? 0 : 10;
    const int myn = yFomFinal.getDataMap().NumMyElements();
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF2::get(numCell, dt, 0.10);
    for (auto i=0; i<myn; i++)
      if (std::abs(yFomFinal[i] - trueY[i+shift]) > 1e-10) checkStr = "FAILED";
  }

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;
}
