
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
  const std::vector<int> corrupt_indices = {1,7,13,19};\
\
  fom_t fomSystem(nFull, corrupt_indices);\
  auto comm = fomSystem.comm();\
  auto full_map = fomSystem.map();\
  fom_state_t fomReferenceState(full_map);\
  pressio::ops::set_zero(fomReferenceState);\
\
  using phi_t = Tpetra::MultiVector<>;\
  phi_t phiFull(full_map, 3);\
  phiFull.getVectorNonConst(0)->putScalar(0.);\
  phiFull.getVectorNonConst(1)->putScalar(1.);\
  phiFull.getVectorNonConst(2)->putScalar(2.);\
\
  auto lv = phiFull.getLocalViewHost();\
  auto mygids = full_map->getMyGlobalIndices();\
  for (std::size_t i=0; i<mygids.extent(0); ++i){\
    if (std::find(corrupt_indices.cbegin(), corrupt_indices.cend(), mygids(i))!=corrupt_indices.cend()){\
      lv(i,0) = -111.;\
      lv(i,1) = -2131.;\
      lv(i,2) =  121.;\
    }\
  }\
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phiFull);\
\
  Eigen::VectorXd romState(3);\
  romState[0]=0.;\
  romState[1]=1.;\
  romState[2]=2.;\


TEST(rom_galerkin_test, const_time_masked_explicit_correctness_tpetra)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t	= TrivialFomOnlyVelocityTpetra;
  MASKED_GALERKIN_COMMON_PART();

  // create masker
  MaskerExplicitTpetra masker(comm, full_map, sample_indices);

  phi_t phiSample(masker.map(), 3);
  phiSample.getVectorNonConst(0)->putScalar(0.);
  phiSample.getVectorNonConst(1)->putScalar(1.);
  phiSample.getVectorNonConst(2)->putScalar(2.);
  ProjectorExplicitTpetra proj(phiSample);

  auto problem = pressio::rom::galerkin::create_masked_explicit_problem(
    pressio::ode::StepScheme::ForwardEuler, 
    fomSystem, decoder, romState, fomReferenceState, proj, masker);

  const scalar_t dt = 1.; 
  const int num_steps = 2;
  ObserverA obs;
  pressio::ode::advance_n_steps_and_observe(problem, romState, 0., dt, num_steps, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 0.);
  EXPECT_DOUBLE_EQ(romState[1], 2611.);
  EXPECT_DOUBLE_EQ(romState[2], 5222.);

  pressio::log::finalize();
}

