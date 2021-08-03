#include <gtest/gtest.h>
#include "pressio_type_traits.hpp"

TEST(type_traits, eigen_dynamic_vector)
{
  using namespace pressio;
  using T = Kokkos::View<double*>;

  using traits = pressio::Traits<T>;
  static_assert(pressio::is_dynamic_vector_kokkos<T>::value, "" );
  static_assert(!pressio::is_static_vector_kokkos<T>::value, "" );
  static_assert(traits::is_static == 0, "" );
  static_assert(std::is_same<typename traits::scalar_type, double>::value, "");
}

TEST(type_traits, eigen_static_vector)
{
  using namespace pressio;
  using T = Kokkos::View<double[5]>;

  using traits = pressio::Traits<T>;
  static_assert(!pressio::is_dynamic_vector_kokkos<T>::value, "" );
  static_assert(pressio::is_static_vector_kokkos<T>::value, "" );
  static_assert(traits::is_static == 1, "" );
  static_assert(std::is_same<typename traits::scalar_type, double>::value, "");
}

