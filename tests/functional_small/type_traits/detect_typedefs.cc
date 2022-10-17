#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

namespace pressio { namespace testing {

#define CHECK2(struct_name, predicate_name, type_name) \
  struct struct_name { \
    using type_name = int; \
  }; \
  static_assert(pressio::predicate_name<struct_name>::value == true, ""); \
  static_assert(pressio::predicate_name<std::ios>::value == false,"");

#define CHECK(struct_name, type_name) CHECK2(struct_name, has_ ## type_name ## def, type_name)

TEST(type_traits, typedef_detect)
{
  CHECK(T1, scalar_type);
  CHECK(T2, state_type);
  CHECK(T3, full_state_type);
  CHECK(T4, reduced_state_type);
  CHECK(T5, basis_matrix_type);
  CHECK(T6, mass_matrix_type);
  CHECK(T7, independent_variable_type);
  CHECK(T8, time_type);
  CHECK(T9, residual_type);
  CHECK(T10, matrix_type);
  CHECK(T11, jacobian_type);
  CHECK(T12, hessian_type);
  CHECK(T13, gradient_type);
  CHECK(T14, discrete_jacobian_type);
  CHECK(T15, discrete_residual_type);
  CHECK(T16, right_hand_side_type);
  CHECK(T17, residual_norm_type);
}

}} // pressio::testing