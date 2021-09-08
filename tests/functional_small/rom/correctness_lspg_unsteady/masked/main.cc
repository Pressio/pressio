
#include <gtest/gtest.h>
#include "../custom_types_specialized_ops.hpp"
#include "../maskers.hpp"
#include "checker.hpp"
#include "foms.hpp"
#include "pressio/rom_lspg.hpp"

template<class T>
void fill_phi(T & phi)
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
  */

  using sc_t = typename pressio::Traits<T>::scalar_type;
  const int nrows = pressio::ops::extent(phi, 0);
  const int ncols = pressio::ops::extent(phi, 1);
  int count = 0;
  for (int i=0; i<nrows; ++i){
    for (int j=0; j<ncols; ++j){
      if (i % 2 == 0){
        phi(i,j) = (sc_t) count++;
      }
      else{
        phi(i,j) = (sc_t) -1;
      }
    }
  }
}

#define COMMON_PART1()\
  using fom_state_t = typename fom_t::state_type;\
  const std::vector<int> rows_to_corrupt_ = {1,3,5,7,9,11,13};\
  const std::vector<int> sample_indices = {0,2,4,6,8,10,12,14};\
  const int nMasked = sample_indices.size();\
  const int N = (int) (rows_to_corrupt_.size() + sample_indices.size());\
  fom_t fomSystem(N, rows_to_corrupt_);\
  fom_state_t fomReferenceState(N);\
  pressio::ops::set_zero(fomReferenceState);\
  Eigen::VectorXd romState(3);\
  romState[0]=0.;\
  romState[1]=1.;\
  romState[2]=2.;\
  phi_t phi(N, 3);\
  fill_phi(phi);\
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

#define COMMON_PART2()\
  const scalar_t dt = 2.;\
  FakeNonLinSolver nonLinSolver(nMasked);\
  ObserverA obs;\
  pressio::ode::advance_n_steps_and_observe(problem.stepper(), romState, 0.,\
              dt, 2, obs, nonLinSolver);\
  std::cout << romState << std::endl;\
  EXPECT_DOUBLE_EQ(romState[0], 4.);\
  EXPECT_DOUBLE_EQ(romState[1], 5.);\
  EXPECT_DOUBLE_EQ(romState[2], 6.);\


TEST(rom_lspg, cont_time_unsteady_masked_correctness_eigen)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomContTimeEigen;
  using scalar_t = typename fom_t::scalar_type;
  using phi_t = Eigen::Matrix<scalar_t, -1,-1>;
  COMMON_PART1();
  MaskerEigen masker(sample_indices);
  auto problem = pressio::rom::lspg::create_masked_unsteady_problem
    (pressio::ode::SteppersE::BDF1, fomSystem, decoder, romState, fomReferenceState, masker);

  COMMON_PART2();

  pressio::log::finalize();
}

TEST(rom_lspg, cont_time_unsteady_masked_correctness_custom_types)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomContTimeCustomTypes;
  using scalar_t = typename fom_t::scalar_type;
  using phi_t = ::pressiotests::MyCustomMatrix<scalar_t>;
  COMMON_PART1();
  MaskerCustomTypes masker(sample_indices);
  auto problem = pressio::rom::lspg::create_masked_unsteady_problem
    (pressio::ode::SteppersE::BDF1, fomSystem, decoder, romState, fomReferenceState, masker);

  COMMON_PART2();

  pressio::log::finalize();
}

TEST(rom_lspg, disc_time_unsteady_masked_correctness_eigen)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomDiscreteTimeEigen;
  using scalar_t = typename fom_t::scalar_type;
  using phi_t = Eigen::Matrix<scalar_t, -1,-1>;

  COMMON_PART1();

  MaskerEigen masker(sample_indices);
  auto problem = pressio::rom::lspg::create_masked_unsteady_problem<2>
    (fomSystem, decoder, romState, fomReferenceState, masker);

  COMMON_PART2();

  pressio::log::finalize();
}


TEST(rom_lspg, disc_time_unsteady_masked_correctness_custom_types)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomDiscreteTimeCustomTypes;
  using scalar_t = typename fom_t::scalar_type;
  using phi_t = ::pressiotests::MyCustomMatrix<scalar_t>;

  COMMON_PART1();

  MaskerCustomTypes masker(sample_indices);
  auto problem = pressio::rom::lspg::create_masked_unsteady_problem<2>
    (fomSystem, decoder, romState, fomReferenceState, masker);

  COMMON_PART2();

  pressio::log::finalize();
}
