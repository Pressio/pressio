
#include <gtest/gtest.h>
#include "./custom_data_types.hpp"
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

// z = beta*z + alpha * A * x
// where x is subscritable as x(i)
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

}}//end namespace pressio::ops

#include "pressio/rom_decoder.hpp"

TEST(rom, linear_decoder_arbitrary_types)
{
  using scalar_t	    = double;
  using fom_state_t	  = pressiotests::MyCustomVector<scalar_t>;
  using decoder_jac_t	= pressiotests::MyCustomMatrix<scalar_t>;

  // create matrix and fill with ones
  decoder_jac_t A(10, 3);
  A.fill(1.);

  // decoder
  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(A);

  // construct reduced state
  // typically, pressio reduced states for ROMs use Eigen or Kokkos (if enabled)
  using rom_state_t = Eigen::VectorXd;
  rom_state_t yRom(A.extent(1));
  yRom.setConstant(2.);

  // apply mapping
  fom_state_t yFom(A.extent(0));
  decoder.applyMapping(yRom, yFom);

  // check solution
  for (std::size_t i=0; i<yFom.extent(0); ++i){
    EXPECT_DOUBLE_EQ(yFom(i), 6.);
  }
}
