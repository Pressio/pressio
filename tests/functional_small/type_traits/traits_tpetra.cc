#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "traits_shared.hpp"

namespace pressio { namespace traits { namespace test {

template <
  typename T,
  int rank,
  typename traits = pressio::Traits<T>
>
void test_tpetra_container()
{
  // traits and shared predicates
  test_container_traits<
    T,
    rank,
    typename T::impl_scalar_type
  >();

  // negative checks (cross-package)
  test_is_not_eigen_container<T>();
  test_is_not_kokkos_container<T>();
  test_is_not_teuchos_container<T>();
}

TEST(tpetra, MVTraits)
{
  using T = Tpetra::MultiVector<
    double,
    int, unsigned int,
    Tpetra::MultiVector<>::node_type
  >;

  // traits and shared predicates
  test_tpetra_container<T, 2>();

  static_assert(pressio::is_multi_vector_tpetra<T>::value,"");

  // negative checks (within Tpetra)
  static_assert(pressio::is_vector_tpetra<T>::value == false, "");
  static_assert(pressio::is_vector_tpetra_block<T>::value == false, "");
  static_assert(pressio::is_multi_vector_tpetra_block<T>::value == false, "");
}

TEST(tpetra, VectorTraits)
{
  using T = Tpetra::Vector<
    double,
    int, unsigned int,
    Tpetra::Vector<>::node_type
  >;

  // traits and shared predicates
  test_tpetra_container<T, 1>();

  // vector predicates
  static_assert(pressio::is_vector_tpetra<T>::value,"");

  // negative checks (within Tpetra)
  static_assert(pressio::is_multi_vector_tpetra<T>::value == false, "");
  static_assert(pressio::is_vector_tpetra_block<T>::value == false, "");
  static_assert(pressio::is_multi_vector_tpetra_block<T>::value == false, "");
}

}}} // pressio::traits::test
