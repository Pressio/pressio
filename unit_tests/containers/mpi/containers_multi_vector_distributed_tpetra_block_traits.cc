#include <gtest/gtest.h>
#include "pressio_containers.hpp"

TEST(containers_multi_vector_distributed_tpetra_block, Traits){
  using namespace pressio;

  using ST = double;
  using LO = int;
  using GO = unsigned long;
  typedef Tpetra::MultiVector<>::node_type NT;
  // typedef Tpetra::Map<LO, GO, NT> map_type;
  using nat_t = Tpetra::Experimental::BlockMultiVector<ST, LO, GO, NT>;
  using mymvec_t = containers::MultiVector<nat_t>;

  static_assert(::pressio::containers::meta::is_multi_vector_tpetra_block<nat_t>::value,"");
  static_assert(::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<mymvec_t>::value,"");

  using mvecTrait = containers::details::traits<mymvec_t>;
  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::scalar_t, ST>();
  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::local_ordinal_t, LO>();

  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::global_ordinal_t,
				unsigned long>();
  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::wrapped_t, nat_t>();
  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::data_map_t,
  				typename nat_t::map_type>();
  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::node_t,
  				typename nat_t::node_type>();
  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::derived_t, mymvec_t>();
  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::communicator_t,
  				Teuchos::RCP<const Teuchos::Comm<int>>
  				>();

  ASSERT_TRUE(mvecTrait::is_vector == false);
  ASSERT_TRUE(mvecTrait::is_multi_vector == true);
  ASSERT_TRUE(mvecTrait::is_shared_mem == false);

  ASSERT_TRUE(mvecTrait::wrapped_multi_vector_identifier
  	      == containers::details::WrappedMultiVectorIdentifier::TpetraBlock);


}//end TEST
