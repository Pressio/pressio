#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

#if defined(USE_EIGEN)
  template<typename Scalar>
  using Type = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

#elif defined(USE_KOKKOS)
  template<typename Scalar>
  using Type = Kokkos::View<Scalar*>;

#else // TODO: use Trilinos types ?
  #define PRESSIO_DISABLE_TEST_SCALARCOMPATIBLE
#endif

using A = Type<double>;
using B = Type<float>;
using C = Type<int>;
using D = Type<bool>;

TEST(traits, all_have_traits_and_same_scalar_2){
  static_assert(  pressio::all_have_traits_and_same_scalar<A, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, B>::value, "");
}

TEST(traits, all_have_traits_and_same_scalar_3){
  static_assert(  pressio::all_have_traits_and_same_scalar<A, A, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, A, B>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, B, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<B, A, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, B, C>::value, "");
}

TEST(traits, all_have_traits_and_same_scalar_4){
  static_assert(  pressio::all_have_traits_and_same_scalar<A, A, A, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, A, A, B>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, A, B, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, B, A, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<B, A, A, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, A, B, B>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, B, A, B>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<B, A, A, B>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, B, B, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<B, A, B, A>::value, "");
  static_assert( !pressio::all_have_traits_and_same_scalar<A, B, C, D>::value, "");
}
