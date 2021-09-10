
#include "foms.hpp"
#include "../checkers.hpp"
#include "pressio/rom_galerkin.hpp"

#define SHARED_LINES()\
  using scalar_t    = typename fom_t::scalar_type;\
  using fom_state_t = typename fom_t::state_type;\
  constexpr int N = 10;\
  fom_t fomSystem(N);\
  fom_state_t fomReferenceState(N);\
  fomReferenceState.setZero();\
  using phi_t = Eigen::MatrixXd;\
  phi_t phi(N, 3);\
  phi.col(0).setConstant(0.);\
  phi.col(1).setConstant(1.);\
  phi.col(2).setConstant(2.);\
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);\
  Eigen::VectorXd romState(3);\
  romState[0]=0.;\
  romState[1]=1.;\
  romState[2]=2.;\


TEST(rom_galerkin, cont_time_default_explicit_correctness_eigen)
{
  /*
    check correctness of Galerkin with Euler forward

    - run two steps: t0 -> t1 -> t2
    - dt = 1.

    - phi in R^{10,3}:
        phi[:,0]=0, phi[:,0]=1, phi[:,0]=2

    - initial romState = [0,1,2]

    - fom velocity f(y,t) always computes y+t

    step0:
      rom state should just be the initial condition

    step1:
      y_fom = phi*[0,1,2]
      rom state = [0,1,2]^T + phi^T f(y_fom, t=0.) = [0,51,102]

    final:
      y_fom = phi*[0,51,102]
      rom state = [0,51,102]^T + phi^T f(y_fom, t=1.) = [0,3601,5202]
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t	= TrivialFomOnlyVelocityEigen;
  SHARED_LINES();

  auto problem = pressio::rom::galerkin::create_default_explicit_problem(
    pressio::ode::SteppersE::ForwardEuler, fomSystem, decoder, romState, fomReferenceState);

  const scalar_t dt = 1.;
  const int num_steps = 2;
  ObserverA obs;
  pressio::ode::advance_n_steps_and_observe(problem.stepper(), romState, 0., dt, num_steps, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 0.);
  EXPECT_DOUBLE_EQ(romState[1], 2611.);
  EXPECT_DOUBLE_EQ(romState[2], 5222.);

  pressio::log::finalize();
}

TEST(rom_galerkin, cont_time_default_implicit_correctness_eigen)
{
  /*
    check correctness of Galerkin with BDF1

    - run two steps: t0 -> t1 -> t2
    - dt = 2.

    - phi in R^{10,3}:
        phi[:,0]=0, phi[:,0]=1, phi[:,0]=2

    - initial romState = [0,1,2]

    - fom velocity f(y,time) always computes y+time
    - fom apply jac (A=J * B) always computes A = B+time
    - fake solver always add 1 to state
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomVelocityAndJacobianEigen;
  SHARED_LINES();

  // using ode_tag = pressio::ode::BDF1;
  auto problem = pressio::rom::galerkin::create_default_implicit_problem(
    pressio::ode::SteppersE::BDF1, fomSystem, decoder, romState, fomReferenceState);
  auto & stepperObj = problem.stepper();

  FakeNonLinSolverContTime nonLinSolver;

  scalar_t dt = 2.;
  pressio::ode::advance_n_steps(stepperObj, romState, 0.0, dt, 2, nonLinSolver);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  pressio::log::finalize();
}


TEST(rom_galerkin, discrete_time_default_implicit_correctness_eigen)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomDiscreteTimeEigen;
  SHARED_LINES();

  auto problem = pressio::rom::galerkin::create_default_problem<2>(fomSystem, decoder, romState, fomReferenceState);
  auto & stepperObj = problem.stepper();

  FakeNonLinSolverForDiscreteTime nonLinSolver;

  scalar_t dt = 2.;
  pressio::ode::advance_n_steps(stepperObj, romState, 0.0, dt, 2, nonLinSolver);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  pressio::log::finalize();
}
