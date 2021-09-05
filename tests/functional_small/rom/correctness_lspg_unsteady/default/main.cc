
#include <gtest/gtest.h>
#include "../custom_types_specialized_ops.hpp"
#include "../checker.hpp"
#include "foms.hpp"
#include "pressio/rom_lspg.hpp"

TEST(rom_lspg, cont_time_unsteady_default_correctness_eigen)
{
  /*
    - phi in R^{10,3}: 
        phi[0,:]=0,1,2
        phi[1,:]=3,4,5
        phi[2,:]=6,7,8
        phi[3,:]=9,10,11
        phi[4,:]=12,13,14
        phi[5,:]=15,16,17
        phi[6,:]=18,19,20
        phi[7,:]=21,22,23

    - initial romState = [0,1,2]

    - fom residual R(y) always computes R[:] = y[:]+1
    - fom applyJac appJac(B) always returns B += 1
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomContTimeEigen;
  using fom_state_t = typename fom_t::state_type;
  using scalar_t = typename fom_t::scalar_type;

  constexpr int N = 8;
  fom_t fomSystem(N);
  fom_state_t fomReferenceState(N);
  fomReferenceState.setZero();

  using phi_t = Eigen::Matrix<scalar_t, -1,-1>;
  phi_t phi(N, 3);
  int count = 0;
  for (int i=0; i<N; ++i){
    for (int j=0; j<3; ++j){
      phi(i,j) = (scalar_t) count++;
    }
  }
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  using tag = pressio::ode::BDF1;
  auto problem = pressio::rom::lspg::create_default_unsteady_problem<tag>(
    fomSystem, decoder, romState, fomReferenceState);

  FakeNonLinSolver nonLinSolver(N);
  const scalar_t dt = 2.;
  ObserverA obs;
  pressio::ode::advance_n_steps_and_observe(problem.stepper(), romState, 0., dt, 1, obs, nonLinSolver);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}

// TEST(rom_lspg, steady_default_correctness_custom_types)
// {
//   // refer to eigen test for description

//   pressio::log::initialize(pressio::logto::terminal);
//   pressio::log::setVerbosity({pressio::log::level::debug});

//   using fom_t	= TrivialFomSteadyCustomTypes;
//   using scalar_t    = typename fom_t::scalar_type;
//   using fom_state_t = typename fom_t::state_type;
//   constexpr int N = 8;
//   fom_t fomSystem(N);
//   fom_state_t fomReferenceState(N);
//   fomReferenceState.fill(0);
//   using phi_t = ::pressiotests::MyCustomMatrix<scalar_t>;
//   phi_t phi(N, 3);
//   int count = 0;
//   for (int i=0; i<N; ++i){
//     for (int j=0; j<3; ++j){
//       phi(i,j) = (scalar_t) count++;
//     }
//   }

//   auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);
//   Eigen::VectorXd romState(3);
//   romState[0]=0.;
//   romState[1]=1.;
//   romState[2]=2.;

//   auto problem = pressio::rom::lspg::create_default_steady_problem(fomSystem, decoder, romState, fomReferenceState);
//   auto & solvableSystem = problem.system();

//   FakeNonLinSolverSteady nonLinSolver(N);
//   nonLinSolver.solve(solvableSystem, romState);
//   std::cout << romState << std::endl;
//   EXPECT_DOUBLE_EQ(romState[0], 2.);
//   EXPECT_DOUBLE_EQ(romState[1], 3.);
//   EXPECT_DOUBLE_EQ(romState[2], 4.);

//   pressio::log::finalize();
// }
