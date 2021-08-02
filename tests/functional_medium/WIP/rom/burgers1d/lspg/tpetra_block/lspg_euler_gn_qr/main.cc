
#include "pressio_rom_lspg.hpp"
#include "pressio_apps.hpp"
#include "utils_tpetra.hpp"

int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dTpetraBlock;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;
  using native_dmat_t = Tpetra::BlockMultiVector<>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;

  std::string checkStr {"PASSED"};

  // scope guard needed for tpetra
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));

    // app object
    constexpr int numCell = 20;
    fom_t appobj( {5.0, 0.02, 0.02}, numCell, Comm);
    scalar_t dt = 0.01;

    // read from file the jacobian of the decoder
    constexpr int romSize = 11;
    // store modes computed before from file
    auto tpw_phi = pressio::rom::test::tpetra::readBasis("basis.txt", romSize, numCell,
         Comm, appobj.getDataMap());
    native_dmat_t tpb_phi(*tpw_phi.data(), *appobj.getDataMap(), 1);
    decoder_jac_t phi(std::move(tpb_phi));
    // create decoder obj
    decoder_t decoderObj(phi);

    // for this problem, my reference state = initial state
    auto & yRef = appobj.getInitialState();

    // define ROM state and initialize to zero (this has to be done)
    lspg_state_t yROM(romSize);
    pressio::ops::fill(yROM, 0.0);

    // define LSPG type
    using ode_tag = pressio::ode::implicitmethods::BDF1;
    auto lspgProblem = pressio::rom::lspg::createDefaultProblemUnsteady<ode_tag>
      (appobj, decoderObj, yROM, yRef);

    // GaussNewton solver
    using rom_jac_t = typename decltype(lspgProblem)::traits::lspg_jacobian_t;
    using qr_solver_type = pressio::qr::QRSolver<rom_jac_t, pressio::qr::TSQR>;
    qr_solver_type qrSolver;

    auto solver = pressio::rom::lspg::createGaussNewtonQRSolver(lspgProblem, yROM, qrSolver);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(4);
    solver.setUpdatingCriterion(pressio::solvers::nonlinear::update::armijo);

    // solve
    pressio::rom::lspg::solveNSequentialMinimizations(lspgProblem, yROM, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.fomStateReconstructorCRef()(yROM);
    auto yFF_v = yFomFinal.data()->getVectorView().getData();

    // this is a reproducing ROM test, so the final reconstructed state
    // has to match the FOM solution obtained with euler, same time-step, for 10 steps
    int shift = (rank==0) ? 0 : 10;
    const int myn = yFomFinal.data()->getMap()->getNodeNumElements();
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
    for (auto i=0; i<myn; i++){
      std::cout << yFF_v[i] << " " << trueY[i+shift] << std::endl;
      if (std::abs(yFF_v[i] - trueY[i+shift]) > 1e-10){
        checkStr = "FAILED";
        break;
      }
    }

  }//tpetra scope

  std::cout << checkStr <<  std::endl;
  return 0;
}
