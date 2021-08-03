#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include "pressio_type_traits.hpp"

TEST(tpetra_block, VectorTraits)
{
  using namespace pressio;

  using ST = double;
  using LO = int;
  using GO = unsigned long;
  typedef Tpetra::BlockVector<>::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  using T = Tpetra::BlockVector<ST, LO, GO, NT>;
  static_assert(pressio::is_vector_tpetra_block<T>::value,"");
  using vecTrait = pressio::Traits<T>;

  ::testing::StaticAssertTypeEq<typename vecTrait::scalar_type, double>();
  ::testing::StaticAssertTypeEq<typename vecTrait::local_ordinal_type,int>();
  ::testing::StaticAssertTypeEq<typename vecTrait::global_ordinal_type,unsigned long>();
  ::testing::StaticAssertTypeEq<typename vecTrait::data_map_type, map_type>();
  ::testing::StaticAssertTypeEq<typename vecTrait::node_type, NT>();
  ::testing::StaticAssertTypeEq<
        typename vecTrait::communicator_type,
				Teuchos::RCP<const Teuchos::Comm<int>>
				>();

  ASSERT_TRUE(vecTrait::rank == 1);
  ASSERT_TRUE(vecTrait::is_shared_mem == false);
  ASSERT_TRUE(vecTrait::is_dynamic == true);
  ASSERT_TRUE(vecTrait::vector_identifier
	      == pressio::VectorIdentifier::TpetraBlock);
}//end TEST

TEST(tpetra_block, MVTraits)
{
  using namespace pressio;

  using ST = double;
  using LO = int;
  using GO = unsigned long;
  typedef Tpetra::MultiVector<>::node_type NT;
  using T = Tpetra::BlockMultiVector<ST, LO, GO, NT>;
  static_assert(::pressio::is_multi_vector_tpetra_block<T>::value,"");

  using mvecTrait = pressio::Traits<T>;
  ::testing::StaticAssertTypeEq<typename mvecTrait::scalar_type, ST>();
  ::testing::StaticAssertTypeEq<typename mvecTrait::local_ordinal_type, LO>();
  ::testing::StaticAssertTypeEq<typename mvecTrait::global_ordinal_type, unsigned long>();
  ::testing::StaticAssertTypeEq<typename mvecTrait::data_map_type,typename T::map_type>();
  ::testing::StaticAssertTypeEq<typename mvecTrait::node_type, typename T::node_type>();
  ::testing::StaticAssertTypeEq<
      typename mvecTrait::communicator_type,
      Teuchos::RCP<const Teuchos::Comm<int>>
      >();

  ASSERT_TRUE(mvecTrait::rank == 2);
  ASSERT_TRUE(mvecTrait::is_shared_mem == false);

  ASSERT_TRUE(mvecTrait::multi_vector_identifier
          == pressio::MultiVectorIdentifier::TpetraBlock);
}//end TEST
