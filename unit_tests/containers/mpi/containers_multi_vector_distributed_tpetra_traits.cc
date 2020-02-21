#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include "pressio_containers.hpp"

TEST(containers_multi_vector_distributed_tpetra, Traits){
  using namespace pressio;

  using ST = double;
  using LO = int;
  using GO = unsigned long;
  typedef Tpetra::MultiVector<>::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  using nat_t = Tpetra::MultiVector<ST, LO, GO, NT>;
  using mymvec_t = containers::MultiVector<nat_t>;

  static_assert(::pressio::containers::meta::is_multi_vector_tpetra<nat_t>::value,"");
  static_assert(::pressio::containers::meta::is_multi_vector_wrapper_tpetra<mymvec_t>::value,"");


  using mvecTrait = containers::details::traits<mymvec_t>;

  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::scalar_t, double>();

  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::local_ordinal_t,int>();

  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::global_ordinal_t,unsigned long>();

  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::wrapped_t, nat_t>();

  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::data_map_t, map_type>();

  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::node_t, NT>();

  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::dot_t, double>();

  ::testing::StaticAssertTypeEq<typename
  				mvecTrait::mag_t, double>();

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
  	      == containers::details::WrappedMultiVectorIdentifier::Tpetra);


}//end TEST
