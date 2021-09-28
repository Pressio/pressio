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

  ::testing::StaticAssertTypeEq<typename
          traits::local_ordinal_type, int>();

  ::testing::StaticAssertTypeEq<typename
          traits::global_ordinal_type, int>();

  ::testing::StaticAssertTypeEq<typename
          traits::data_map_type, Epetra_BlockMap>();

  ::testing::StaticAssertTypeEq<typename
          traits::size_type, int>();

  ::testing::StaticAssertTypeEq<typename
          traits::communicator_type, Epetra_Comm>();
}

TEST(epetra, VectorTraits)
{
  using T = Epetra_Vector;
  static_assert(pressio::is_vector_epetra<T>::value,"");
  static_assert(pressio::Traits<T>::vector_identifier == pressio::VectorIdentifier::Epetra);
  test_epetra_container<T, 1>();
}

TEST(eped_epetra, MVTraits)
{
  using T = Epetra_MultiVector;
  static_assert(pressio::is_multi_vector_epetra<T>::value,"");
  static_assert(pressio::Traits<T>::multi_vector_identifier == pressio::MultiVectorIdentifier::Epetra);
  test_epetra_container<T, 2>();
}