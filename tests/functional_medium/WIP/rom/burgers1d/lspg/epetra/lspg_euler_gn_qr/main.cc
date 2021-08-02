
#include "pressio_rom_lspg.hpp"
#include "pressio_apps.hpp"
#include "utils_epetra.hpp"

int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

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
  scalar_t dt = 0.01;

  // read from file the jacobian of the decoder
  constexpr int romSize = 11;
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::rom::test::epetra::readBasis("basis.txt", romSize, numCell,
					 Comm, appobj.getDataMap());
  const int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;
  // create decoder obj
  decoder_t decoderObj(phi);

  // for this problem, my reference state = initial state
  auto & yRef = appobj.getInitialState();

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  pressio::ops::fill( yROM, 0.0);

  // define LSPG type
  using ode_tag  = pressio::ode::implicitmethods::BDF1;
  // using lspg_problem = typename pressio::rom::lspg::composeDefaultProblem<
  //   ode_tag, fom_t, decoder_t, lspg_state_t>::type;
  // lspg_problem lspgProblem(appobj, decoderObj, yROM, yRef);
  auto lspgProblem = pressio::rom::lspg::createDefaultProblemUnsteady<ode_tag>
    (appobj, decoderObj, yROM, yRef);

  // GaussNewton solver
  using rom_jac_t = typename decltype(lspgProblem)::traits::lspg_jacobian_t;
  using qr_solver_type = pressio::qr::QRSolver<rom_jac_t, pressio::qr::TSQR>;
  qr_solver_type qrSolver;

  auto solver = pressio::rom::lspg::createGaussNewtonQRSolver(lspgProblem, yROM, qrSolver);
  solver.setUpdatingCriterion(pressio::solvers::nonlinear::update::armijo);
  solver.setTolerance(1e-13);
  // I know this should converge in few iters at every step
  solver.setMaxIterations(5);

  // solve
  pressio::rom::lspg::solveNSequentialMinimizations(lspgProblem, yROM, 0.0, dt, 10, solver);

  // compute the fom corresponding to our rom final state
  auto yFomFinal = lspgProblem.fomStateReconstructorCRef()(yROM);
  yFomFinal.data()->Print(std::cout << std::setprecision(14));

  // this is a reproducing ROM test, so the final reconstructed state
  // has to match the FOM solution obtained with euler, same time-step, for 10 steps
  {
    int shift = (rank==0) ? 0 : 10;
    const int myn = yFomFinal.data()->Map().NumMyElements();
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
    for (auto i=0; i<myn; i++)
      if (std::abs(yFomFinal(i) - trueY[i+shift]) > 1e-10) checkStr = "FAILED";
  }

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;
}
