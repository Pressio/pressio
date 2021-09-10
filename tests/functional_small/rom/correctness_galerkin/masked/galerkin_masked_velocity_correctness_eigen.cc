
#include "../../custom_data_types.hpp"
#include "foms.hpp"
#include "projectors.hpp"
#include "maskers.hpp"
#include "../checkers.hpp"
#include "pressio/rom_galerkin.hpp"

#define MASKED_GALERKIN_COMMON_PART() \
  using scalar_t    = typename fom_t::scalar_type;\
  using fom_state_t = typename fom_t::state_type;\
\
  constexpr int nFull = 20;\
  const std::vector<int> sample_indices = {0,2,4,6,8,10,12,14,16,18};\
  const int nMasked = sample_indices.size();\
  /* corrupt indices are those that we mess up on purpose */\
  const std::vector<int> corrupt_indices = {1,7,13,19};\
\
  fom_t fomSystem(nFull, corrupt_indices);\
  fom_state_t fomReferenceState(nFull);\
  fomReferenceState.setZero();\
\
  /* the decoder must be valid to reconstruct fhe full fom */\
  using phi_t = Eigen::MatrixXd;\
  phi_t phiFull(nFull, 3);\
  phiFull.col(0).setConstant(0.);\
  phiFull.col(1).setConstant(1.);\
  phiFull.col(2).setConstant(2.);\
  /* corrupting the entries in phiFull is only way we can ensure the masking works properly */ \
  phiFull.row(1).setConstant(-111.);\
  phiFull.row(7).setConstant(111.);\
  phiFull.row(11).setConstant(423.);\
  phiFull.row(17).setConstant(-21.);\
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phiFull);\
\
  Eigen::VectorXd romState(3);\
  romState[0]=0.;\
  romState[1]=1.;\
  romState[2]=2.;\
  /* projector must be applicable to the *masked* operand*/\
  /* so we need to use only certain rows of phi*/\
  phi_t phiSample(nMasked, 3);\
  for (int i = 0; i < nMasked; ++i){\
    phiSample(i, 0) = phiFull(sample_indices[i],0);\
    phiSample(i, 1) = phiFull(sample_indices[i],1);\
    phiSample(i, 2) = phiFull(sample_indices[i],2);\
  }\


TEST(rom_galerkin_test, const_time_masked_explicit_correctness_eigen)
{
  /*
    check correctness of masked Galerkin 
    this test, *after the mask is applied*, 
    is doing the same thing the default galerkin does. 
    so that is why the values are the same as the default test.
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomOnlyVelocityEigen;
  MASKED_GALERKIN_COMMON_PART();

  MaskerExplicitEigen masker(sample_indices);
  ProjectorExplicitEigen proj(phiSample);

  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_masked_explicit_problem(
    pressio::ode::SteppersE::ForwardEuler, fomSystem, decoder, romState, fomReferenceState, proj, masker);

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

TEST(rom_galerkin_test, const_time_masked_implicit_correctness_eigen)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomVelocityAndJacobianEigen;
  MASKED_GALERKIN_COMMON_PART();

  MaskerImplicitEigen masker(sample_indices);
  ProjectorImplicitEigen proj(phiSample);

  auto problem = pressio::rom::galerkin::create_masked_implicit_problem(
    pressio::ode::SteppersE::BDF1, fomSystem, decoder, romState, fomReferenceState, proj, masker);
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

TEST(rom_galerkin_test, discrete_time_masked_implicit_correctness_eigen)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomDiscreteTimeEigen;
  MASKED_GALERKIN_COMMON_PART();

  MaskerImplicitEigen masker(sample_indices);
  ProjectorImplicitEigen proj(phiSample);

  auto problem = pressio::rom::galerkin::create_masked_problem<2>(
    fomSystem, decoder, romState, fomReferenceState, proj, masker);
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


