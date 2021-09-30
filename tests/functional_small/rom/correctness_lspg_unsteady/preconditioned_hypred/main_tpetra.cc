
#include <gtest/gtest.h>
#include "../custom_types_specialized_ops.hpp"
#include "../preconditioners.hpp"
#include "checker.hpp"
#include "updater.hpp"
#include "foms.hpp"
#include "pressio/rom_lspg.hpp"

TEST(rom_lspg, cont_time_unsteady_hyperreduced_correctness_tpetra)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t	= TrivialFomVelocityAndJacobianTpetra;
  using scalar_t = typename fom_t::scalar_type;
  using fom_state_t = typename fom_t::state_type;

  int N = 20;
  fom_t fomSystem(N);
  auto comm = fomSystem.comm();
  int rank = comm->getRank();
  auto map = fomSystem.map();
  fom_state_t fomReferenceState(map);
  pressio::ops::set_zero(fomReferenceState);

  using phi_t = Tpetra::MultiVector<>;
  phi_t phi(map, 3);
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  PreconditionerTpetra prec(rank);
  auto problem = pressio::rom::lspg::create_prec_hyperreduced_unsteady_problem
    (pressio::ode::StepScheme::BDF1, fomSystem, decoder, romState, fomReferenceState, prec);

  const scalar_t dt = 2.;
  FakeNonLinSolverTpetra nonLinSolver;
  ObserverA obs;
  pressio::ode::advance_n_steps_and_observe(problem.stepper(), romState, 0.,
					    dt, 2, obs, nonLinSolver);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  pressio::log::finalize();
}
