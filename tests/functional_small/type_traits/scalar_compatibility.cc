#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

#if defined(PRESSIO_ENABLE_TPL_EIGEN) // use Eigen types for testing
  template<typename Scalar>
  using Type = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
#elif defined(PRESSIO_ENABLE_TPL_KOKKOS) // use Kokkos types for testing
  template<typename Scalar>
  using Type = Kokkos::View<Scalar*>;
#else // TODO: use Trilinos types ?
  #define PRESSIO_DISABLE_TEST_SCALARCOMPATIBLE
#endif

#ifndef PRESSIO_DISABLE_TEST_SCALARCOMPATIBLE

using A = Type<double>;
using B = Type<float>;
using C = Type<int>;
using D = Type<bool>;

TEST(containers_meta, two_vector_scalar_compatible){
  static_assert(  pressio::are_scalar_compatible<A, A>::value, "");
  static_assert( !pressio::are_scalar_compatible<A, B>::value, "");
}

TEST(containers_meta, three_vector_scalar_compatible){
  static_assert(  pressio::are_scalar_compatible<A, A, A>::value, "");
  static_assert( !pressio::are_scalar_compatible<A, A, B>::value, "");
  static_assert( !pressio::are_scalar_compatible<A, B, A>::value, "");
  static_assert( !pressio::are_scalar_compatible<B, A, A>::value, "");
  static_assert( !pressio::are_scalar_compatible<A, B, C>::value, "");
}

TEST(containers_meta, four_vector_scalar_compatible){
  static_assert(  pressio::are_scalar_compatible<A, A, A, A>::value, "");

  static_assert( !pressio::are_scalar_compatible<A, A, A, B>::value, "");
  static_assert( !pressio::are_scalar_compatible<A, A, B, A>::value, "");
  static_assert( !pressio::are_scalar_compatible<A, B, A, A>::value, "");
  static_assert( !pressio::are_scalar_compatible<B, A, A, A>::value, "");

  static_assert( !pressio::are_scalar_compatible<A, A, B, B>::value, "");
  static_assert( !pressio::are_scalar_compatible<A, B, A, B>::value, "");
  static_assert( !pressio::are_scalar_compatible<B, A, A, B>::value, "");
  static_assert( !pressio::are_scalar_compatible<A, B, B, A>::value, "");
  static_assert( !pressio::are_scalar_compatible<B, A, B, A>::value, "");

  static_assert( !pressio::are_scalar_compatible<A, B, C, D>::value, "");
}

#endif