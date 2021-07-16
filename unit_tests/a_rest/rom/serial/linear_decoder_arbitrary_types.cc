
#include <gtest/gtest.h>
#include "pressio_rom.hpp"

template <typename scalar_t>
struct MyOps
{
  using z_t = std::vector<scalar_t>;
  using A_t = std::vector<std::vector<scalar_t>>;

  // z = beta*z + alpha * A * x
  // where x is something that is subscritable as x(i)
  template< typename x_t>
  void product(pressio::nontranspose,
	       scalar_t alpha,
	       const A_t & A,
	       const x_t & x,
	       scalar_t beta,
	       z_t & z) const
  {
    // obviously not efficient, just for demonstration
    for (std::size_t i=0; i<A.size(); ++i)
    {
      z[i] += beta*z[i];
      for (std::size_t j=0; j<A[i].size(); ++j){
	z[i] += alpha*A[i][j]*x(j);
      }
    }
  }
};

TEST(rom_decoder, linearArbitraryType)
{
  // *** define some types ***
  // here we assume your FOM application uses an Eigen vector for the state
  // and an Eigen matrix as the type for the Jacobian
  using scalar_t	   = double;
  using native_fom_state_t = std::vector<scalar_t>;
  using native_phi_t	   = std::vector<std::vector<scalar_t>>;

  // the wrapped types
  // what happens in pressio: std::vector is treated as unknwon type by pressio
  // so effectively pressio::containers::Vector is labeled as an "arbitrary" type
  using fom_state_t	= pressio::containers::Vector<native_fom_state_t>;
  using decoder_jac_t	= pressio::containers::DenseMatrix<native_phi_t>;

  // *** fill phi ***
  // create a native phi with 10 rows and 3 columns and fill with ones
  native_phi_t phiNative(6);
  for (auto & iRow : phiNative){
    iRow.resize(2, 1.);
  }

  // *** construct decoder  ***
  using ops_t = MyOps<scalar_t>;
  using decoder_t = pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t, ops_t>;
  // Need to pass the native matrix and an object that knows
  // how to compute the operations (see MyOps at the top)
  ops_t ops;
  decoder_t decoder(phiNative, ops);

  // *** construct reduced state  ***
  // typically, pressio reduced states for ROMs use Eigen or Kokkos (if enabled)
  using rom_state_t = pressio::containers::Vector<Eigen::VectorXd>;
  rom_state_t yRom(2);
  // set yRom = 2.
  pressio::ops::fill(yRom, 2.);

  // *** apply mapping ***
  fom_state_t yFom(6);
  decoder.applyMapping(yRom, yFom);

  // *** check solution ***
  // yFom should be = [4. 4. .... 4.]
  for (auto i=0; i<6; ++i) EXPECT_DOUBLE_EQ(yFom(i), 4.);
}
