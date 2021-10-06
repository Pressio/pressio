#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "traits_shared.hpp"

#define CHECK_TRAIT2(TRAIT, TYPE) \
::testing::StaticAssertTypeEq<typename traits::TRAIT, typename TYPE>();
#define CHECK_TRAIT(TRAIT) CHECK_TRAIT2(TRAIT, T::TRAIT)

template <
  typename T,
  typename traits = pressio::Traits<T>
>
void test_tpetra_extra_traits()
{
  CHECK_TRAIT(dual_view_type);
  CHECK_TRAIT(dot_type);
  CHECK_TRAIT(mag_type);
}

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

  CHECK_TRAIT(local_ordinal_type);
  CHECK_TRAIT(global_ordinal_type);
  CHECK_TRAIT2(data_map_type, T::map_type);
  CHECK_TRAIT2(communicator_type, pressio::impl::TrilinosCommType<T>::type);
  CHECK_TRAIT(node_type);
  CHECK_TRAIT(device_type);
  CHECK_TRAIT2(device_mem_space_type, T::device_type::memory_space);
  CHECK_TRAIT2(device_exec_space_type, T::device_type::execution_space);
  CHECK_TRAIT2(host_mem_space_type, Kokkos::HostSpace::memory_space);
  CHECK_TRAIT2(host_exec_space_type, Kokkos::HostSpace::execution_space);
}

TEST(tpetra, MVTraits)
{
  using T = Tpetra::MultiVector<
    double,
    int, unsigned int,
    Tpetra::MultiVector<>::node_type
  >;
  static_assert(pressio::is_multi_vector_tpetra<T>::value,"");
  ASSERT_TRUE(pressio::Traits<T>::multi_vector_identifier
      == pressio::MultiVectorIdentifier::Tpetra);
  test_tpetra_container<T, 2>();
  test_tpetra_extra_traits<T>();
}

TEST(tpetra, VectorTraits)
{
  using T = Tpetra::Vector<
    double,
    int, unsigned int,
    Tpetra::Vector<>::node_type
  >;
  static_assert(pressio::is_vector_tpetra<T>::value,"");
  ASSERT_TRUE(pressio::Traits<T>::vector_identifier
        == pressio::VectorIdentifier::Tpetra);
  test_tpetra_container<T, 1>();
  test_tpetra_extra_traits<T>();
}

TEST(tpetra_block, VectorTraits)
{
  using T = Tpetra::BlockVector<
    double,
    int, unsigned int,
    Tpetra::BlockVector<>::node_type
  >;
  static_assert(::pressio::is_vector_tpetra_block<T>::value,"");
  ASSERT_TRUE(pressio::Traits<T>::vector_identifier
        == pressio::VectorIdentifier::TpetraBlock);
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
  ASSERT_TRUE(pressio::Traits<T>::multi_vector_identifier
        == pressio::MultiVectorIdentifier::TpetraBlock);
  test_tpetra_container<T, 2>();
}
