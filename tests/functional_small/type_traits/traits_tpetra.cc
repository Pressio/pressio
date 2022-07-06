#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "traits_shared.hpp"

#define CHECK_TRAIT2(TRAIT, TYPE) \
::testing::StaticAssertTypeEq<typename traits::TRAIT, typename TYPE>();
#define CHECK_TRAIT(TRAIT) CHECK_TRAIT2(TRAIT, T::TRAIT)


template <
  typename T,
  int rank,
  typename traits = pressio::Traits<T>
>
void test_tpetra_container()
{
  test_container_traits<
    T,
    pressio::PackageIdentifier::Trilinos,
    rank,
    false, /* shared mem */
    true, /* dynamic */
    typename T::impl_scalar_type,
    typename T::local_ordinal_type,
    typename T::global_ordinal_type
  >();
}

TEST(tpetra, MVTraits)
{
  using T = Tpetra::MultiVector<
    double,
    int, unsigned int,
    Tpetra::MultiVector<>::node_type
  >;
  static_assert(pressio::is_multi_vector_tpetra<T>::value,"");
  test_tpetra_container<T, 2>();
}

TEST(tpetra, VectorTraits)
{
  using T = Tpetra::Vector<
    double,
    int, unsigned int,
    Tpetra::Vector<>::node_type
  >;
  static_assert(pressio::is_vector_tpetra<T>::value,"");
  test_tpetra_container<T, 1>();
}

TEST(tpetra_block, VectorTraits)
{
  using T = Tpetra::BlockVector<
    double,
    int, unsigned int,
    Tpetra::BlockVector<>::node_type
  >;
  static_assert(::pressio::is_vector_tpetra_block<T>::value,"");
  test_tpetra_container<T, 1>();
}

TEST(tpetra_block, MVTraits)
{
  using T = Tpetra::BlockMultiVector<
    double,
    int, unsigned int,
    Tpetra::MultiVector<>::node_type
  >;
  static_assert(pressio::is_multi_vector_tpetra_block<T>::value,"");
  test_tpetra_container<T, 2>();
}
