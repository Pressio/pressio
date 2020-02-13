#include <gtest/gtest.h>
// #include <Tpetra_Map.hpp>
#include "CONTAINERS_ALL"

TEST(containers_sparse_matrix_tpetra, Traits){
  using namespace pressio;

  using ST = double;
  using LO = Tpetra::CrsMatrix<>::local_ordinal_type;
  using GO = Tpetra::CrsMatrix<>::global_ordinal_type;
  using NT = Tpetra::CrsMatrix<>::node_type;
  using map_type =  Tpetra::Map<LO, GO, NT>;

  using nat_t = Tpetra::CrsMatrix<>;
  STATIC_ASSERT_IS_MATRIX_SPARSE_DISTRIBUTED_TPETRA(nat_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(nat_t);
  STATIC_ASSERT_IS_NOT_MULTIVECTOR_EPETRA(nat_t);
  STATIC_ASSERT_IS_NOT_VECTOR_TPETRA(nat_t);

  using mymat_t = containers::Matrix<nat_t>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(mymat_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(mymat_t);
  STATIC_ASSERT_IS_NOT_MULTIVECTOR_EPETRA(mymat_t);
  STATIC_ASSERT_IS_NOT_MULTIVECTOR_TPETRA(mymat_t);
  STATIC_ASSERT_IS_CONTAINERS_MATRIX_WRAPPER(mymat_t);

  static_assert(containers::meta::is_sparse_matrix_wrapper_tpetra<mymat_t>::value, " ");

  using matTrait = containers::details::traits<mymat_t>;

  ::testing::StaticAssertTypeEq<typename
  				matTrait::scalar_t, ST>();

  ::testing::StaticAssertTypeEq<typename
  				matTrait::local_ordinal_t,LO>();

  ::testing::StaticAssertTypeEq<typename
			matTrait::global_ordinal_t,GO>();

  ::testing::StaticAssertTypeEq<typename
  				matTrait::wrapped_t, nat_t>();

  ::testing::StaticAssertTypeEq<typename
  				matTrait::row_map_t, map_type>();
  ::testing::StaticAssertTypeEq<typename
  				matTrait::col_map_t, map_type>();
  ::testing::StaticAssertTypeEq<typename
  				matTrait::range_map_t, map_type>();
  ::testing::StaticAssertTypeEq<typename
  				matTrait::domain_map_t, map_type>();

  ::testing::StaticAssertTypeEq<typename
  				matTrait::node_t, NT>();
  ::testing::StaticAssertTypeEq<typename
  				matTrait::mag_t, double>();
  ::testing::StaticAssertTypeEq<typename
  				matTrait::derived_t, mymat_t>();

  ::testing::StaticAssertTypeEq<typename
  				matTrait::communicator_t,
  				Teuchos::RCP<const Teuchos::Comm<int>>
  				>();

  ASSERT_TRUE(matTrait::is_matrix == true);
  ASSERT_TRUE(matTrait::is_vector == false);
  ASSERT_TRUE(matTrait::is_multi_vector == false);
  ASSERT_TRUE(matTrait::is_shared_mem == false);

  ASSERT_TRUE(matTrait::wrapped_matrix_identifier
  	      == containers::details::WrappedMatrixIdentifier::SparseTpetra);


}//end TEST
