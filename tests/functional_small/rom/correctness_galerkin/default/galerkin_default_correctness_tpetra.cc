
#include "foms.hpp"
#include "../checkers.hpp"
#include "pressio/rom_galerkin.hpp"

TEST(rom_galerkin, cont_time_default_explicit_correctness_tpetra)
{
  // refer to eigen test for description

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t	= TrivialFomOnlyVelocityTpetra;
  using scalar_t    = typename fom_t::scalar_type;
  using fom_state_t	= typename fom_t::state_type;

  constexpr int N = 10;
  fom_t fomSystem(N);
  auto comm = fomSystem.comm();
  auto map = fomSystem.map();
  fom_state_t fomReferenceState(map);
  pressio::ops::set_zero(fomReferenceState);

  using phi_t = Tpetra::MultiVector<>;
  phi_t phi(map, 3);
  phi.getVectorNonConst(0)->putScalar(0.);
  phi.getVectorNonConst(1)->putScalar(1.);
  phi.getVectorNonConst(2)->putScalar(2.);
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  using ode_tag = pressio::ode::ForwardEuler;
  auto problem = pressio::rom::galerkin::create_default_problem<ode_tag>(fomSystem, decoder, romState, fomReferenceState);

  const scalar_t dt = 1.; 
  const int num_steps = 2;
  ObserverA obs;
  pressio::ode::advance_n_steps_and_observe(problem.stepper(), romState, 0., dt, num_steps, obs);

  EXPECT_DOUBLE_EQ(romState[0], 0.);
  EXPECT_DOUBLE_EQ(romState[1], 2611.);
  EXPECT_DOUBLE_EQ(romState[2], 5222.);

  pressio::log::finalize();
}


TEST(rom_galerkin, cont_time_default_implicit_correctness_tpetra)
{
  // refer to eigen test for description

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomVelocityAndJacobianTpetra;
  using scalar_t    = typename fom_t::scalar_type;
  using fom_state_t = typename fom_t::state_type;

  constexpr int N = 10;
  fom_t fomSystem(N);
  auto comm = fomSystem.comm();
  auto map = fomSystem.map();
  fom_state_t fomReferenceState(map);
  pressio::ops::set_zero(fomReferenceState);

  using phi_t = Tpetra::MultiVector<>;
  phi_t phi(map, 3);
  phi.getVectorNonConst(0)->putScalar(0.);
  phi.getVectorNonConst(1)->putScalar(1.);
  phi.getVectorNonConst(2)->putScalar(2.);
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  using ode_tag = pressio::ode::BDF1;
  auto problem = pressio::rom::galerkin::create_default_problem<ode_tag>(fomSystem, decoder, romState, fomReferenceState);
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
