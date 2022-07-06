#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "traits_shared.hpp"

template <
  typename T,
  int rank,
  typename traits = pressio::Traits<T>
>
void test_epetra_container()
{
  test_container_traits<
    T,
    pressio::PackageIdentifier::Trilinos,
    rank,
    false, /* shared mem */
    true, /* dynamic */
    double, /* scalar */
    int /* ordinal */
  >();
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