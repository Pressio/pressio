#include <gtest/gtest.h>
#include "CORE_VECTOR"
#include "CORE_MULTI_VECTOR"

TEST(core_multi_vector_distributed_tpetra_block, Traits){
  using namespace rompp;

  using ST = double;
  using LO = int;
  using GO = unsigned long;
  typedef Tpetra::MultiVector<>::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  using nat_t = Tpetra::Experimental::BlockMultiVector<ST, LO, GO, NT>;

  STATIC_ASSERT_IS_MULTIVECTOR_TPETRA_BLOCK(nat_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(nat_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(nat_t);
  STATIC_ASSERT_IS_NOT_MULTIVECTOR_EPETRA(nat_t);
  STATIC_ASSERT_IS_NOT_VECTOR_TPETRA(nat_t);

  using mymvec_t = core::MultiVector<nat_t>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(mymvec_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(mymvec_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(mymvec_t);
  STATIC_ASSERT_IS_NOT_MULTIVECTOR_EPETRA(mymvec_t);
  STATIC_ASSERT_IS_NOT_MULTIVECTOR_TPETRA(mymvec_t);
  STATIC_ASSERT_IS_CORE_MULTI_VECTOR_WRAPPER(mymvec_t);

  using mvecTrait = core::details::traits<mymvec_t>;
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
  	      == core::details::WrappedMultiVectorIdentifier::TpetraBlock);


}//end TEST
