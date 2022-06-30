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

TEST(traits, all_have_traits)
{
  static_assert(  pressio::all_have_traits<A>::value, "");
  static_assert(  pressio::all_have_traits<B>::value, "");
  static_assert(  pressio::all_have_traits<C>::value, "");
  static_assert(  pressio::all_have_traits<D>::value, "");

  static_assert(  !pressio::all_have_traits<double>::value, "");
  static_assert(  !pressio::all_have_traits<int>::value, "");

  static_assert(  pressio::all_have_traits<A,B>::value, "");
  static_assert(  !pressio::all_have_traits<A,double>::value, "");
  static_assert(  pressio::all_have_traits<A,B,C>::value, "");
  static_assert(  pressio::all_have_traits<A,B,C,D>::value, "");
  static_assert(  !pressio::all_have_traits<A,int,B,C,D>::value, "");
}
