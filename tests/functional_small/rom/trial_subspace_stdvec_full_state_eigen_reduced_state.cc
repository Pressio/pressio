
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

namespace pressio{
template<>
struct Traits<std::vector<double>>{
  using scalar_type = double;
  static constexpr int rank = 1;
};
}

#include "pressio/ops.hpp"

namespace pressio{ namespace ops{
template<class T>
std::vector<T> clone(const std::vector<T> & o){
  return o;
}

template<class T>
void fill(std::vector<T> & o, T value){
  std::for_each(o.begin(), o.end(), [=](T & entry){ entry = value; });
}

// y = beta*y + alpha*A*x
void product(::pressio::nontranspose,
	     double alpha, const Eigen::MatrixXd & A,
	     const Eigen::VectorXd & x,
	     double beta,
	     std::vector<double> & y)
{
  for (std::size_t i=0; i<y.size(); ++i){
    double sum = 0;
    for (int j=0; j<A.cols(); ++j){
      sum += A(i,j) * x[j];
    }
    y[i] = beta*y[i] + alpha*sum;
  }
}

void update(std::vector<double> & y,
	    double a,
	    const std::vector<double> & x,
	    double b)
{
  for (std::size_t i=0; i<y.size(); ++i){
    y[i] = a*y[i] + b*x[i];
  }
}

}} // end namespace pressio::ops

#include "pressio/rom_subspaces.hpp"

TEST(rom, trial_subspace_construct_1)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi;
  using reduced_state_type = Eigen::VectorXd;
  std::vector<double> shift(15);
  auto space = create_trial_column_subspace<reduced_state_type>(phi, shift, false);
  (void) space;
}

TEST(rom, trial_subspace_construct_2)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  using reduced_state_type = Eigen::VectorXd;
  std::vector<double> shift(15);
  auto space = create_trial_column_subspace<reduced_state_type>(basis_t(), shift, false);
}

TEST(rom, trial_subspace_construct_3)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi;
  using reduced_state_type = Eigen::VectorXd;
  std::vector<double> shift(15);
  auto space = create_trial_column_subspace<reduced_state_type>(phi, std::move(shift), false);
  (void) space;
}

TEST(rom, trial_subspace_construct_4)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  using reduced_state_type = Eigen::VectorXd;
  std::vector<double> shift(15);
  auto space = create_trial_column_subspace<reduced_state_type>(basis_t(), std::move(shift), false);
  (void) space;
}

TEST(rom, trial_subspace_create_reduced_state)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi(15,3);

  using reduced_state_type = Eigen::VectorXd;
  std::vector<double> shift(15);
  auto space = create_trial_column_subspace<reduced_state_type>(std::move(phi), shift, false);
  auto a = space.createReducedState();
  EXPECT_TRUE(a.size() == 3);
  for (int i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], 0.);
  }
}

TEST(rom, trial_subspace_create_full_state)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi(15,3);

  using reduced_state_type = Eigen::VectorXd;
  std::vector<double> shift(15);
  auto space = create_trial_column_subspace<reduced_state_type>(std::move(phi), shift, false);
  auto a = space.createFullState();
  EXPECT_TRUE(a.size() == 15);

  for (decltype(a.size()) i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], 0.);
  }
}

TEST(rom, trial_subspace_map_from_reduced_state)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi = basis_t::Random(15,3);

  using reduced_state_type = Eigen::VectorXd;
  std::vector<double> shift(15);
  auto space = create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  Eigen::VectorXd gold = phi * latState;

  auto a = space.createFullState();
  space.mapFromReducedState(latState, a);
  for (decltype(a.size()) i=0; i<a.size(); ++i){
    EXPECT_NEAR(a[i], gold[i], 1e-15);
  }
}

TEST(rom, trial_subspace_create_full_from_reduced)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi = basis_t::Random(15,3);

  using reduced_state_type = Eigen::VectorXd;
  std::vector<double> shift(15);
  auto space = create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  Eigen::VectorXd gold = phi * latState;

  auto a = space.createFullStateFromReducedState(latState);
  for (decltype(a.size()) i=0; i<a.size(); ++i){
    EXPECT_NEAR(a[i], gold[i], 1e-15);
  }
}

TEST(rom, affine_trial_subspace_construct_1)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);
  pressio::ops::fill(shift, 0.);
  auto space = create_trial_column_subspace<reduced_state_type>(basis_t(), shift, true);
}

TEST(rom, affine_trial_subspace_construct_2)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);
  pressio::ops::fill(shift, 0.);
  auto space = create_trial_column_subspace<reduced_state_type>(phi, shift, true);
  (void) space;
}

TEST(rom, affine_trial_subspace_create_reduced_state)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi(15,3);

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);
  pressio::ops::fill(shift, 0.);
  auto space = create_trial_column_subspace<reduced_state_type>(std::move(phi), shift, true);
  auto a = space.createReducedState();
  EXPECT_TRUE(a.size() == 3);
  for (decltype(a.size()) i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], 0.);
  }
}

TEST(rom, affine_trial_subspace_create_full_state)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi(15,3);

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);
  pressio::ops::fill(shift, 0.);

  auto space = create_trial_column_subspace<reduced_state_type>(std::move(phi), shift, true);
  auto a = space.createFullState();
  EXPECT_TRUE(a.size() == 15);

  for (decltype(a.size()) i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], 0.);
  }
}

TEST(rom, affine_trial_subspace_map_from_reduced_state)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi = basis_t::Random(15,3);

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);
  pressio::ops::fill(shift, 2.);

  auto space = create_trial_column_subspace<reduced_state_type>(phi, shift, true);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  Eigen::VectorXd gold = phi * latState;

  auto a = space.createFullState();
  space.mapFromReducedState(latState, a);
  for (decltype(a.size()) i=0; i<a.size(); ++i){
    EXPECT_NEAR(a[i], gold[i] + 2., 1e-15);
  }
}

TEST(rom, affine_trial_subspace_create_full_from_reduced)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi = basis_t::Random(15,3);

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);
  pressio::ops::fill(shift, 2.);

  auto space = create_trial_column_subspace<reduced_state_type>(phi, shift, true);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  Eigen::VectorXd gold = phi * latState;

  auto a = space.createFullStateFromReducedState(latState);
  for (decltype(a.size()) i=0; i<a.size(); ++i){
    EXPECT_NEAR(a[i], gold[i] + 2., 1e-15);
  }
}

TEST(rom, affine_trial_subspace_view_basis)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi = basis_t::Random(15,3);

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);

  auto space = create_trial_column_subspace<reduced_state_type>(phi, shift, true);

  reduced_state_type latState1 = reduced_state_type::Random(3);
  auto & J = space.basis();
  EXPECT_TRUE( J.isApprox(phi) );
}


TEST(rom, trial_subspace_shift_is_zero_if_nonaffine)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi;
  using reduced_state_type = Eigen::VectorXd;
  const double fillvalueshift = 1145.4;
  std::vector<double> shift(15, fillvalueshift);
  auto space = create_trial_column_subspace<reduced_state_type>(phi, shift, false);
  const auto & shiftStored = space.translationVector();
  ASSERT_TRUE( std::all_of(shiftStored.cbegin(), shiftStored.cend(), [](auto v){ return v==0; }) );
}
