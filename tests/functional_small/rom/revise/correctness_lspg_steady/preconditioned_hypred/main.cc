

#include <gtest/gtest.h>
#include "../custom_types_specialized_ops.hpp"
#include "../checker.hpp"
#include "foms.hpp"
#include "../preconditioners.hpp"
#include "pressio/rom_lspg.hpp"

TEST(rom_lspg, steady_preconitioned_hypred_correctness_eigen)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t	= TrivialFomSteadyEigen;
  using fom_state_t	= typename fom_t::state_type;
  using scalar_t = typename fom_t::scalar_type;

  const int nstencil = 15;
  const int nSample  = 8;
  fom_t fomSystem(nSample);
  fom_state_t fomReferenceState(nstencil);
  fomReferenceState.setZero();

  using phi_t = Eigen::Matrix<scalar_t, -1,-1>;
  phi_t phi(nstencil, 3);
  int count = 0;
  for (int i=0; i<nstencil; ++i){
    for (int j=0; j<3; ++j){
      if (i % 2 == 0){
        phi(i,j) = (scalar_t) count++;
      }
      else{
        phi(i,j) = (scalar_t) -1;
      }
    }
  }
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  PreconditionerSteadyEigen prec;

  auto problem = pressio::rom::lspg::create_hyperreduced_steady_problem(fomSystem, decoder, romState, fomReferenceState, prec);

  FakeNonLinSolverSteady nonLinSolver(nSample);
  nonLinSolver.solve(problem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}

TEST(rom_lspg, steady_preconitioned_hypred_correctness_custom_types)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomSteadyCustomTypes;
  using fom_state_t = typename fom_t::state_type;
  using scalar_t = typename fom_t::scalar_type;

  const int nstencil = 15;
  const int nSample  = 8;
  fom_t fomSystem(nSample);
  fom_state_t fomReferenceState(nstencil);
  fomReferenceState.fill(0);

  using phi_t = ::pressiotests::MyCustomMatrix<scalar_t>;
  phi_t phi(nstencil, 3);
  int count = 0;
  for (int i=0; i<nstencil; ++i){
    for (int j=0; j<3; ++j){
      if (i % 2 == 0){
        phi(i,j) = (scalar_t) count++;
      }
      else{
        phi(i,j) = (scalar_t) -1;
      }
    }
  }
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  PreconditionerSteadyCustomTypes<scalar_t> prec;

  auto problem = pressio::rom::lspg::create_hyperreduced_steady_problem(fomSystem, decoder, romState, fomReferenceState, prec);

  FakeNonLinSolverSteady nonLinSolver(nSample);
  nonLinSolver.solve(problem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}
