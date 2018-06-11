
#include <gtest/gtest.h>
#include "vector/core_vector_serial_eigen.hpp"
#include "vector/core_vector_serial_stdlib.hpp"
#include "vector/core_vector_serial_userdefined.hpp"
#include "vector/core_vector_meta.hpp"
//#include "core_static_assert_definitions.hpp"

TEST(coreVector, EigenVector)
{
  using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  static_assert( std::is_same<eigvec_t, Eigen::VectorXd>::value, "--");

  // create an eigen vector
  eigvec_t e_v1;

  // create wrapper for eigen vector
  using myvec_t = core::vector<eigvec_t>;
  myvec_t m_v1;
  myvec_t m_v1;
  myvec_t m_v1;
  
}
