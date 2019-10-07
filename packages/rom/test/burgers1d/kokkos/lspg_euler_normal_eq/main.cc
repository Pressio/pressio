
#include "UTILS_ALL"
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_kokkos.hpp"

int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dKokkos;
  using scalar_t	= typename fom_t::scalar_type;

  // using exe_space = typename fom_t::execution_space;
  using native_state_t_d = typename fom_t::state_type_d;
  using native_state_t_h = typename fom_t::state_type_h;
  using native_mv_t_d = typename fom_t::mv_d;
  // using native_mv_t_h = typename fom_t::mv_h;

  // device lspg state type
  using lspg_state_d_t	= pressio::containers::Vector<native_state_t_d>;

  // device decoder jacobian type
  using decoder_jac_d_t	= pressio::containers::MultiVector<native_mv_t_d>;
  // host decoder jacobian type
  // using decoder_jac_h_t	= pressio::containers::MultiVector<native_mv_t_h>;

  // device decoder type
  using decoder_d_t	= pressio::rom::LinearDecoder<decoder_jac_d_t>;

  std::string checkStr {"PASSED"};
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  Kokkos::initialize (argc, argv);
  {
    // app object
    constexpr int numCell = 20;
    fom_t appobj({{5.0, 0.02, 0.02}}, numCell);
    appobj.setup();
    const scalar_t dt = 0.01;
    constexpr auto t0 = zero;

    constexpr int romSize = 11;

    // create/read jacobian of the decoder
    decoder_jac_d_t phi("phi", numCell, 11);
    pressio::rom::test::kokkos::readBasis("basis.txt", romSize, numCell, *phi.data());
    if( phi.numVectors() != romSize ) return 0;

    // create decoder obj
    decoder_d_t decoderObj(phi);

    // for this problem, my reference state = initial state
    auto yRef = appobj.getInitialState();

    // define ROM state
    lspg_state_d_t yROM("yRom", romSize);

    // define LSPG type
    constexpr auto ode_case  = pressio::ode::ImplicitEnum::Euler;
    using lspg_problem_t = pressio::rom::DefaultLSPGTypeGenerator<
      fom_t, ode_case, decoder_d_t, lspg_state_d_t>;
    pressio::rom::LSPGUnsteadyProblemGenerator<lspg_problem_t> lspgProblem(appobj,
    									   yRef,
    									   decoderObj,
    									   yROM, t0);
    using lspg_stepper_t = typename lspg_problem_t::lspg_stepper_t;

    // linear solver
    using hessian_t  = pressio::containers::Matrix<typename fom_t::mv_d>;
    using solver_tag   = pressio::solvers::linear::direct::getrs;
    using linear_solver_t = pressio::solvers::direct::KokkosDirect<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using gnsolver_t   = pressio::solvers::iterative::GaussNewton<lspg_stepper_t,
    								  linear_solver_t>;
    gnsolver_t solver(lspgProblem.getStepperRef(), yROM, linSolverObj);
    solver.setTolerance(1e-14);
    // I know this should converge in few iters at every step
    solver.setMaxIterations(2);

    // integrate in time
    pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    const auto yFomFinal_d = lspgProblem.getFomStateReconstructorCRef()(yROM);

    // create a host mirror for yFomFinal
    native_state_t_h yFomFinal_h("yFF_h", numCell);
    Kokkos::deep_copy(yFomFinal_h, *yFomFinal_d.data());

    // this is a reproducing ROM test, so the final reconstructed state
    // has to match the FOM solution obtained with euler, same time-step, for 10 steps
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
    for (auto i=0; i<numCell; i++){
      std::cout << std::setprecision(15) << yFomFinal_h(i) << " " << trueY[i] << std::endl;
      if (std::abs(yFomFinal_h(i) - trueY[i]) > 1e-10) checkStr = "FAILED";
    }

    auto n1 = ::pressio::containers::ops::norm2(yFomFinal_d);
    std::cout << n1 << std::endl;

    std::cout << checkStr <<  std::endl;
  }
  Kokkos::finalize();
  return 0;
}
