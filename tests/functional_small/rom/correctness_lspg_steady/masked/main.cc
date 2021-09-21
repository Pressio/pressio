
#include <gtest/gtest.h>
#include "../custom_types_specialized_ops.hpp"
#include "../checker.hpp"
#include "foms.hpp"
#include "../maskers.hpp"
#include "pressio/rom_lspg.hpp"

TEST(rom_lspg, steady_masked_correctness_eigen)
{
  /*
    - phi in R^{10,3}: 
        phi[0,:]=0,1,2
        phi[1,:]=-1,-1,-1
        phi[2,:]=3,4,5
        phi[3,:]=-1,-1,-1
        phi[4,:]=6,7,8
        phi[5,:]=-1,-1,-1
        phi[6,:]=9,10,11
        phi[7,:]=-1,-1,-1
        phi[8,:]=12,13,14
        phi[9,:]=-1,-1,-1
        phi[10,:]=15,16,17
        phi[11,:]=-1,-1,-1
        phi[12,:]=18,19,20
        phi[13,:]=-1,-1,-1
        phi[14,:]=21,22,23

    - initial romState = [0,1,2]

    - fom residual R(y) always computes R[:] = y[:]+1
    - fom applyJac appJac(B) always returns B += 1
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomSteadyEigen;
  using fom_state_t = typename fom_t::state_type;
  using scalar_t = typename fom_t::scalar_type;

  const std::vector<int> indices_to_corrupt_ = {1,3,5,7,9,11,13};

  constexpr int N = 15;
  fom_t fomSystem(N, indices_to_corrupt_);
  fom_state_t fomReferenceState(N);
  fomReferenceState.setZero();

  using phi_t = Eigen::Matrix<scalar_t, -1,-1>;
  phi_t phi(N, 3);
  int count = 0;
  for (int i=0; i<N; ++i){
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

  const std::vector<int> sample_indices = {0,2,4,6,8,10,12,14};
  const int nMasked = sample_indices.size();
  MaskerSteadyEigen masker(sample_indices);

  auto problem = pressio::rom::lspg::create_masked_steady_problem(fomSystem, decoder, romState, fomReferenceState, masker);
  auto & solvableSystem = problem.system();

  FakeNonLinSolverSteady nonLinSolver(nMasked);
  nonLinSolver.solve(solvableSystem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}

TEST(rom_lspg, steady_masked_correctness_custom_types)
{

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomSteadyCustomTypes;
  using fom_state_t = typename fom_t::state_type;
  using scalar_t = typename fom_t::scalar_type;

  const std::vector<int> indices_to_corrupt_ = {1,3,5,7,9,11,13};

  constexpr int N = 15;
  fom_t fomSystem(N, indices_to_corrupt_);
  fom_state_t fomReferenceState(N);
  fomReferenceState.fill(0);

  using phi_t = ::pressiotests::MyCustomMatrix<scalar_t>;
  phi_t phi(N, 3);
  int count = 0;
  for (int i=0; i<N; ++i){
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

  const std::vector<int> sample_indices = {0,2,4,6,8,10,12,14};
  const int nMasked = sample_indices.size();
  MaskerSteadyCustomTypes masker(sample_indices);

  auto problem = pressio::rom::lspg::create_masked_steady_problem(fomSystem, decoder, romState, fomReferenceState, masker);
  auto & solvableSystem = problem.system();

  FakeNonLinSolverSteady nonLinSolver(nMasked);
  nonLinSolver.solve(solvableSystem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}
