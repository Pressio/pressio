
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/rom_subspaces.hpp"
#include "fixtures/tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorGlobSize15Fixture, trial_subspace_construct_1)
{
  auto phi = ::pressio::ops::clone(*myMv_);

  using namespace pressio::rom;
  using basis_t = mvec_t;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = vec_t;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(std::move(phi));
  (void) space;
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, trial_subspace_construct_2)
{
  auto phi = ::pressio::ops::clone(*myMv_);

  using namespace pressio::rom;
  using basis_t = mvec_t;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = vec_t;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(phi);
  (void) space;
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, trial_subspace_create_reduced_state)
{
  using namespace pressio::rom;
  using basis_t = mvec_t;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = vec_t;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(*myMv_);

  auto a = space.createReducedState();
  EXPECT_TRUE(a.size() == 4);
  for (int i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], 0.);
  }
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, trial_subspace_create_full_state)
{
  using namespace pressio::rom;
  using basis_t = mvec_t;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = vec_t;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(*myMv_);

  auto a = space.createFullState();
  EXPECT_TRUE( ::pressio::ops::extent(a,0) == 15 );

  auto a_h = a.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
  for (int i=0; i<a.getLocalLength(); ++i){
    EXPECT_DOUBLE_EQ(a_h(i,0), 0.);
  }
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, trial_subspace_map_from_reduced_state)
{

  using namespace pressio::rom;
  using basis_t = mvec_t;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = vec_t;
  auto phi = *myMv_;
  phi.randomize();
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(phi);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  auto phi_h = myMv_->getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
  Eigen::VectorXd gold(phi_h.extent(0));
  for (int i=0; i<phi_h.extent(0); ++i){
    gold[i] = 0.;
    for (int j=0; j<phi_h.extent(1); ++j){
      gold[i] += phi_h(i,j) * latState(j);
    }
  }

  auto a = space.createFullState();
  space.mapFromReducedState(latState, a);
  auto a_h = a.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
  for (int i=0; i<a.getLocalLength(); ++i){
    EXPECT_NEAR(a_h(i,0), gold[i], 1e-12);
  }
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, trial_subspace_create_full_from_reduced)
{

  using namespace pressio::rom;
  using basis_t = mvec_t;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = vec_t;
  auto phi = *myMv_;
  phi.randomize();
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(phi);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  auto phi_h = myMv_->getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
  Eigen::VectorXd gold(phi_h.extent(0));
  for (int i=0; i<phi_h.extent(0); ++i){
    gold[i] = 0.;
    for (int j=0; j<phi_h.extent(1); ++j){
      gold[i] += phi_h(i,j) * latState(j);
    }
  }

  auto a = space.createFullStateFromReducedState(latState);
  auto a_h = a.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
  for (int i=0; i<a.getLocalLength(); ++i){
    EXPECT_NEAR(a_h(i,0), gold[i], 1e-12);
  }
}
