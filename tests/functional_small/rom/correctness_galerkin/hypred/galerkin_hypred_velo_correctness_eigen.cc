
#include <gtest/gtest.h>
#include "../../custom_data_types.hpp"
#include "foms.hpp"
#include "projectors.hpp"
#include "../checkers.hpp"
#include "pressio/rom_galerkin.hpp"

#define HYPRED_VELO_GALERKIN_COMMON_PART() \
  using scalar_t    = typename fom_t::scalar_type;\
  using fom_state_t = typename fom_t::state_type;\
\
  const int nstencil = 20;\
  const int nSample  = 10;\
  fom_t fomSystem(nSample);\
  fom_state_t fomReferenceState(nstencil);\
  fomReferenceState.setZero();\
\
  using phi_t = Eigen::MatrixXd;\
  phi_t phi(nstencil, 3);\
  phi.col(0).setConstant(0.);\
  phi.col(1).setConstant(1.);\
  phi.col(2).setConstant(2.);\
  phi.row(0).setConstant(-111.);\
  phi.row(2).setConstant(-111.);\
  phi.row(4).setConstant(111.);\
  phi.row(6).setConstant(423.);\
  phi.row(8).setConstant(-21.);\
  phi.row(10).setConstant(423.);\
  phi.row(12).setConstant(-21.);\
  phi.row(14).setConstant(423.);\
  phi.row(16).setConstant(-21.);\
  phi.row(18).setConstant(-21.);\
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);\
\
  Eigen::VectorXd romState(3);\
  romState[0]=0.;\
  romState[1]=1.;\
  romState[2]=2.;\
  /* projector must be applicable to the *sample* operand */\
  phi_t matForProj(nSample, 3);\
  matForProj.col(0).setConstant(0.);\
  matForProj.col(1).setConstant(1.);\
  matForProj.col(2).setConstant(2.);\
  ProjectorExplicitEigen proj(matForProj);\


TEST(rom_galerkin, cont_time_hypred_explicit_correctness_eigen)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomOnlyVelocityEigen;
  HYPRED_VELO_GALERKIN_COMMON_PART();

  using ode_tag = pressio::ode::ForwardEuler;
  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_hyperreduced_problem<ode_tag>(
    fomSystem, decoder, romState, fomReferenceState, proj);

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

TEST(rom_galerkin, cont_time_hypred_implicit_correctness_eigen)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomVelocityAndJacobianEigen;
  HYPRED_VELO_GALERKIN_COMMON_PART();

  using ode_tag = pressio::ode::BDF1;
  auto problem = pressio::rom::galerkin::create_hyperreduced_problem<ode_tag>(
    fomSystem, decoder, romState, fomReferenceState, proj);
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

TEST(rom_galerkin, discrete_time_hypred_implicit_correctness_eigen)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomDiscreteTimeEigen;
  HYPRED_VELO_GALERKIN_COMMON_PART();

  auto problem = pressio::rom::galerkin::create_hyperreduced_problem<2>(
    fomSystem, decoder, romState, fomReferenceState, proj);
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
