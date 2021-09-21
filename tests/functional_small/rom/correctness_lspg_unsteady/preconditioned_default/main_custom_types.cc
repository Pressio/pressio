
#include <gtest/gtest.h>
#include "../custom_types_specialized_ops.hpp"
#include "../preconditioners.hpp"
#include "checker.hpp"
#include "foms.hpp"
#include "pressio/rom_lspg.hpp"

template<class T>
void fill_phi(T & phi)
{
  using sc_t = typename pressio::Traits<T>::scalar_type;
  const int nrows = pressio::ops::extent(phi, 0);
  const int ncols = pressio::ops::extent(phi, 1);
  int count = 0;
  for (int i=0; i<nrows; ++i){
    for (int j=0; j<ncols; ++j){
      phi(i,j) = (sc_t) count++;
    }
  }
}

TEST(rom_lspg, cont_time_unsteady_default_correctness_custom_types)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomContTimeCustomTypes;
  using fom_state_t = typename fom_t::state_type;
  using scalar_t = typename fom_t::scalar_type;

  constexpr int N = 8;
  fom_t fomSystem(N);
  fom_state_t fomReferenceState(N);
  fomReferenceState.fill(0);

  using phi_t = ::pressiotests::MyCustomMatrix<scalar_t>;
  phi_t phi(N, 3);
  fill_phi(phi);
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  PreconditionerCustomTypes<scalar_t> prec;
  auto problem = pressio::rom::lspg::create_default_unsteady_problem
    (pressio::ode::StepScheme::BDF1, fomSystem, decoder, romState, fomReferenceState, prec);

  const scalar_t dt = 2.;
  FakeNonLinSolver<phi_t> nonLinSolver(N, phi, dt);
  ObserverA obs;
  pressio::ode::advance_n_steps_and_observe(problem.stepper(), romState, 0.,
					    dt, 2, obs, nonLinSolver);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  pressio::log::finalize();
}

TEST(rom_lspg, disc_time_prec_default_correctness_custom_types)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomDiscreteTimeCustomTypes;
  using fom_state_t = typename fom_t::state_type;
  using scalar_t = typename fom_t::scalar_type;

  constexpr int N = 8;
  fom_t fomSystem(N);
  fom_state_t fomReferenceState(N);
  fomReferenceState.fill(0);

  using phi_t = ::pressiotests::MyCustomMatrix<scalar_t>;
  phi_t phi(N, 3);
  fill_phi(phi);
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  PreconditionerCustomTypes<scalar_t> prec;
  auto problem = pressio::rom::lspg::create_default_unsteady_problem<2>
    (fomSystem, decoder, romState, fomReferenceState, prec);

  const scalar_t dt = 2.;
  FakeNonLinSolver<phi_t> nonLinSolver(N, phi, dt);
  ObserverA obs;
  pressio::ode::advance_n_steps_and_observe(problem.stepper(), romState, 0.,
					    dt, 2, obs, nonLinSolver);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  pressio::log::finalize();
}
