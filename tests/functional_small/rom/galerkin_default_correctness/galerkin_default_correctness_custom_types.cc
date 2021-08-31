
#include <gtest/gtest.h>
#include "common.hpp"
#include "pressio/ops.hpp"

namespace pressio{ 

template<class ScalarType> 
struct Traits<pressiotests::MyCustomVector<ScalarType>>{
  using scalar_type = ScalarType;
};

template<class ScalarType> 
struct Traits<pressiotests::MyCustomMatrix<ScalarType>>{
  using scalar_type = ScalarType;
};

namespace ops{

template<class ScalarType> 
std::size_t extent(pressiotests::MyCustomVector<ScalarType> & object, int i){
    return object.extent(i);
}

template<class ScalarType> 
std::size_t extent(pressiotests::MyCustomMatrix<ScalarType> & object, int i){
    return object.extent(i);
}

template<class ScalarType> 
void set_zero(pressiotests::MyCustomVector<ScalarType> & object){
  object.fill(0);
}

template<class ScalarType> 
void set_zero(pressiotests::MyCustomMatrix<ScalarType> & object){
  object.fill(0);
}

template<class ScalarType> 
void deep_copy(pressiotests::MyCustomVector<ScalarType> & dest, 
               const pressiotests::MyCustomVector<ScalarType> & src){
  dest = src;
}

template<class ScalarType> 
void deep_copy(pressiotests::MyCustomMatrix<ScalarType> & dest, 
               const pressiotests::MyCustomMatrix<ScalarType> & src){
  dest = src;
}

template<class ScalarType> 
pressiotests::MyCustomVector<ScalarType> clone(const pressiotests::MyCustomVector<ScalarType> & src){
  return pressiotests::MyCustomVector<ScalarType>(src.extent(0));
}

template<class ScalarType> 
pressiotests::MyCustomMatrix<ScalarType> clone(const pressiotests::MyCustomMatrix<ScalarType> & src){
  return pressiotests::MyCustomMatrix<ScalarType>(src.extent(0), src.extent(1));
}

template<class ScalarType> 
void update(pressiotests::MyCustomVector<ScalarType> & v,        const ScalarType a,
            const pressiotests::MyCustomVector<ScalarType> & v1, const ScalarType b)
{
  for (std::size_t i=0; i< v.extent(0); ++i){
    v(i) = a*v(i) + b*v1(i);
  }
}

// z = beta*z + alpha * A * x
// where x is indexable as x(i)
template<class x_t, class ScalarType>
void product(pressio::nontranspose,
             ScalarType alpha,
             const pressiotests::MyCustomMatrix<ScalarType> & A,
             const x_t & x,
             ScalarType beta,
             pressiotests::MyCustomVector<ScalarType> & z)
{
  // obviously not efficient, just for demonstration
  for (std::size_t i=0; i<A.extent(0); ++i)
  {
    z(i) = beta*z(i);
    for (std::size_t j=0; j<A.extent(1); ++j){
      z(i) += alpha*A(i,j)*x(j);
    }
  }
}

// z = beta*z + alpha * A^T * x
// where z is indexable as z(i)
template<class z_t, class ScalarType>
void product(pressio::transpose,
             ScalarType alpha,
             const pressiotests::MyCustomMatrix<ScalarType> & A,
             const pressiotests::MyCustomVector<ScalarType> & x,
             ScalarType beta,
             z_t & z)
{
  // obviously not efficient, just for demonstration
  for (std::size_t k=0; k<A.extent(1); ++k)
  {
    z(k) = beta*z(k);
    for (std::size_t i=0; i<A.extent(0); ++i){
      z(k) += alpha*A(i,k)*x(i);
    }
  }
}

// C = beta*C + alpha * A^T * B
// where C is indexable as C(i,j)
template<class C_t, class ScalarType>
void product(pressio::transpose,
             pressio::nontranspose,
             ScalarType alpha,
             const pressiotests::MyCustomMatrix<ScalarType> & A,
             const pressiotests::MyCustomMatrix<ScalarType> & B,
             ScalarType beta,
             C_t & C)
{
  for (std::size_t i=0; i<A.extent(1); ++i){
    for (std::size_t j=0; j<B.extent(1); ++j)
    {
      C(i,j) *= beta;
      for (std::size_t k=0; k<A.extent(0); ++k){
        C(i,j) += alpha*A(k,i)*B(k,j);
      }
    }
  }
}

}}//end namespace pressio::ops

#include "pressio/rom_galerkin.hpp"

TEST(rom_galerkin, default_explicit_correctness_custom_types)
{
  // refer to eigen test for description

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t	= TrivialFomOnlyVelocityCustomTypes;
  using scalar_t    = typename fom_t::scalar_type;
  using fom_state_t	= typename fom_t::state_type;

  constexpr int N = 10;
  fom_t fomSystem(N);
  fom_state_t fomReferenceState(N);
  fomReferenceState.fill(0);

  using phi_t = ::pressiotests::MyCustomMatrix<scalar_t>;
  phi_t phi(N, 3);
  for (std::size_t i=0; i<phi.extent(0); ++i){
    phi(i,0) = 0;
    phi(i,1) = 1;
    phi(i,2) = 2;
  }
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  using ode_tag = pressio::ode::ForwardEuler;
  auto problem = pressio::rom::galerkin::create_default_problem<ode_tag>(fomSystem, decoder, romState, fomReferenceState);

  const scalar_t dt = 1.; 
  const int num_steps = 2;
  ObserverA obs;
  pressio::ode::advance_n_steps_and_observe(problem.stepper(), romState, 0., dt, num_steps, obs);

  EXPECT_DOUBLE_EQ(romState[0], 0.);
  EXPECT_DOUBLE_EQ(romState[1], 2611.);
  EXPECT_DOUBLE_EQ(romState[2], 5222.);

  pressio::log::finalize();
}

TEST(rom_galerkin, default_implicit_correctness_custom_types)
{
  // refer to eigen test for description

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using fom_t = TrivialFomVelocityAndJacobianCustomTypes;
  using scalar_t    = typename fom_t::scalar_type;
  using fom_state_t = typename fom_t::state_type;

  constexpr int N = 10;
  fom_t fomSystem(N);
  fom_state_t fomReferenceState(N);
  fomReferenceState.fill(0);

  using phi_t = ::pressiotests::MyCustomMatrix<scalar_t>;
  phi_t phi(N, 3);
  for (std::size_t i=0; i<phi.extent(0); ++i){
    phi(i,0) = 0;
    phi(i,1) = 1;
    phi(i,2) = 2;
  }
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  Eigen::VectorXd romState(3);
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  using ode_tag = pressio::ode::BDF1;
  auto problem = pressio::rom::galerkin::create_default_problem<ode_tag>(fomSystem, decoder, romState, fomReferenceState);
  auto & stepperObj = problem.stepper();

  FakeNonLinSolver nonLinSolver;

  scalar_t dt = 2.;
  pressio::ode::advance_n_steps(stepperObj, romState, 0.0, dt, 2, nonLinSolver);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  pressio::log::finalize();
}
