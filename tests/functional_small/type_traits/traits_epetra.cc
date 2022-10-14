#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "traits_shared.hpp"

namespace pressio { namespace traits { namespace test {

template <
  typename T,
  int rank,
  typename traits = pressio::Traits<T>
>
void test_epetra_container()
{
  test_container_traits<
    T,
    rank,
    double /* scalar */
  >();

  // negative
  test_is_not_eigen_container<T>();
  test_is_not_teuchos_container<T>();
  test_is_not_tpetra_container<T>();
  test_is_not_kokkos_container<T>();
}

TEST(epetra, VectorTraits)
{
  using T = Epetra_Vector;
  static_assert(pressio::is_vector_epetra<T>::value,"");
  test_epetra_container<T, 1>();
}

TEST(eped_epetra, MVTraits)
{
  using T = Epetra_MultiVector;
  static_assert(pressio::is_multi_vector_epetra<T>::value,"");
  test_epetra_container<T, 2>();
}

}}} // pressio::traits::test