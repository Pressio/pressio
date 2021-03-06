#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include "pressio_containers.hpp"

TEST(containers_vector_distributed_tpetra_block, Traits){
  using namespace pressio;

  using ST = double;
  using LO = int;
  using GO = unsigned long;
  typedef Tpetra::BlockVector<>::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  using natV_t = Tpetra::BlockVector<ST, LO, GO, NT>;
  static_assert(::pressio::containers::predicates::is_vector_tpetra_block<natV_t>::value,"");

  using myvec_t = containers::Vector<natV_t>;
  static_assert(::pressio::containers::predicates::is_vector_wrapper_tpetra_block<myvec_t>::value,"");

  using vecTrait = containers::details::traits<myvec_t>;

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::scalar_t, double>();
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::local_ordinal_t,int>();
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::global_ordinal_t,unsigned long>();
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::data_map_t, map_type>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::wrapped_t, natV_t>();
  ::testing::StaticAssertTypeEq<typename
				vecTrait::node_t, NT>();
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::communicator_t,
				Teuchos::RCP<const Teuchos::Comm<int>>
				>();

  ASSERT_TRUE(vecTrait::rank == 1);
  ASSERT_TRUE(vecTrait::is_shared_mem == false);
  ASSERT_TRUE(vecTrait::is_dynamic == true);
  ASSERT_TRUE(vecTrait::wrapped_vector_identifier
	      == containers::details::WrappedVectorIdentifier::TpetraBlock);


}//end TEST
