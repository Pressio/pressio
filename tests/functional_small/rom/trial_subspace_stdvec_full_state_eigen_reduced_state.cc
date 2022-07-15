
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

namespace pressio{
template<>
struct Traits<std::vector<double>>{
  using scalar_type = double;
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
  for (int i=0; i<y.size(); ++i){
    double sum = 0;
    for (int j=0; j<A.cols(); ++j){
      sum += A(i,j) * x[j];
    }
    y[i] = beta*y[i] + alpha*sum;
  }
}

// y = alpha*A*x
template<
  class ResultType,
  mpl::enable_if_t< std::is_same<ResultType, std::vector<double>>::value, int > = 0
  >
ResultType product(::pressio::nontranspose,
		   double alpha, const Eigen::MatrixXd & A,
		   const Eigen::VectorXd & x)
{
  std::vector<double> y(A.rows());
  for (int i=0; i<y.size(); ++i){
    double sum = 0;
    for (int j=0; j<A.cols(); ++j){
      sum += A(i,j) * x[j];
    }
    y[i] = alpha*sum;
  }
  return y;
}

void update(std::vector<double> & y,
	    double a,
	    const std::vector<double> & x,
	    double b)
{
  for (int i=0; i<y.size(); ++i){
    y[i] = a*y[i] + b*x[i];
  }
}

}} // end namespace pressio::ops

#include "pressio/rom_subspaces.hpp"

TEST(rom, trial_subspace_construct_1)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(basis_t());
}

TEST(rom, trial_subspace_construct_2)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(phi);
}

TEST(rom, trial_subspace_create_reduced_state)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi(15,3);

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(std::move(phi));
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
  using full_state_type = std::vector<double>;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(std::move(phi));
  auto a = space.createFullState();
  EXPECT_TRUE(a.size() == 15);

  for (int i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], 0.);
  }
}

TEST(rom, trial_subspace_map_from_reduced_state)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi = basis_t::Random(15,3);

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(phi);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  Eigen::VectorXd gold = phi * latState;

  auto a = space.createFullState();
  space.mapFromReducedState(latState, a);
  for (int i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], gold[i]);
  }
}

TEST(rom, trial_subspace_create_full_from_reduced)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi = basis_t::Random(15,3);

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  auto space = create_trial_subspace<reduced_state_type, full_state_type>(phi);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  Eigen::VectorXd gold = phi * latState;

  auto a = space.createFullStateFromReducedState(latState);
  for (int i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], gold[i]);
  }
}

TEST(rom, affine_trial_subspace_construct_1)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);
  auto space = create_affine_trial_subspace<reduced_state_type>(basis_t(), shift);
}

TEST(rom, affine_trial_subspace_construct_2)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi;
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);
  auto space = create_affine_trial_subspace<reduced_state_type>(phi, shift);
}

TEST(rom, affine_trial_subspace_create_reduced_state)
{
  using namespace pressio::rom;

  using basis_t = Eigen::MatrixXd;
  basis_t phi(15,3);

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = std::vector<double>;
  full_state_type shift(15);
  auto space = create_affine_trial_subspace<reduced_state_type>(std::move(phi), shift);
  auto a = space.createReducedState();
  EXPECT_TRUE(a.size() == 3);
  for (int i=0; i<a.size(); ++i){
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

  auto space = create_affine_trial_subspace<reduced_state_type>(std::move(phi), shift);
  auto a = space.createFullState();
  EXPECT_TRUE(a.size() == 15);

  for (int i=0; i<a.size(); ++i){
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

  auto space = create_affine_trial_subspace<reduced_state_type>(phi, shift);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  Eigen::VectorXd gold = phi * latState;

  auto a = space.createFullState();
  space.mapFromReducedState(latState, a);
  for (int i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], gold[i] + 2.);
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

  auto space = create_affine_trial_subspace<reduced_state_type>(phi, shift);

  auto latState = space.createReducedState();
  latState = reduced_state_type::Random(latState.size());
  Eigen::VectorXd gold = phi * latState;

  auto a = space.createFullStateFromReducedState(latState);
  for (int i=0; i<a.size(); ++i){
    EXPECT_DOUBLE_EQ(a[i], gold[i] + 2.);
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

  auto space = create_affine_trial_subspace<reduced_state_type>(phi, shift);

  reduced_state_type latState1 = reduced_state_type::Random(3);
  auto & J = space.viewBasis();
  EXPECT_TRUE( J.isApprox(phi) );
}
