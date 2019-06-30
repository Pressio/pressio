
#include "epetra_only_fixtures.hpp"
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h"
#include "CONTAINERS_ALL"


TEST_F(epetraSparseMatR7MultiVectorR9C4Fixture,
       CrsMatrixMVecEpetraProduct){

  using mymat_w_t = rompp::containers::Matrix<Epetra_CrsMatrix>;
  using mymv_w_t = rompp::containers::MultiVector<Epetra_MultiVector>;  
  ASSERT_TRUE(rompp::containers::meta::is_sparse_matrix_wrapper_epetra<mymat_w_t>::value);
  ASSERT_TRUE(rompp::containers::meta::is_multi_vector_wrapper_epetra<mymv_w_t>::value);

  // fill data
  fillCrsMatrix();
  fillMultiVector();
  
  // we need to call fillcomplete considering rectangular matrix because
  // the multivector has 9 rows, while matrix has 7 rows
  // C = A * mv so:
  // use as range map the row map of the matrix 
  const auto & rangeMap = *dataMapSM_; // the range map = map of sparse matrix
  const auto & domainMap = *dataMapMV_; // the domain map is map of multivector
  sm_->fillingIsCompleted(domainMap, rangeMap);
  //  sm_->data()->Print(std::cout);

  // compute prod
  auto C = rompp::containers::ops::product(*sm_, *mv_);
  C.data()->Print(std::cout);
  
  if (rank_==0){
    ASSERT_EQ(C.localLength(),3);
    EXPECT_DOUBLE_EQ( C(0,0), 0. );
    EXPECT_DOUBLE_EQ( C(0,1), 6. );
    EXPECT_DOUBLE_EQ( C(0,2), 1. );
    EXPECT_DOUBLE_EQ( C(0,3), 0. );

    EXPECT_DOUBLE_EQ( C(1,0), 0. );
    EXPECT_DOUBLE_EQ( C(1,1), 5. );
    EXPECT_DOUBLE_EQ( C(1,2), 2. );
    EXPECT_DOUBLE_EQ( C(1,3), 0. );

    EXPECT_DOUBLE_EQ( C(2,0), 2. );
    EXPECT_DOUBLE_EQ( C(2,1), 5. );
    EXPECT_DOUBLE_EQ( C(2,2), 1. );
    EXPECT_DOUBLE_EQ( C(2,3), 0. );    
  }

  if (rank_==1){
    ASSERT_EQ(C.localLength(),2);
    EXPECT_DOUBLE_EQ( C(0,0), 6. );
    EXPECT_DOUBLE_EQ( C(0,1), 12. );
    EXPECT_DOUBLE_EQ( C(0,2), 6. );
    EXPECT_DOUBLE_EQ( C(0,3), 0. );

    EXPECT_DOUBLE_EQ( C(1,0), 1. );
    EXPECT_DOUBLE_EQ( C(1,1), 2. );
    EXPECT_DOUBLE_EQ( C(1,2), 0. );
    EXPECT_DOUBLE_EQ( C(1,3), 0. );
  }

  if (rank_==2){
    ASSERT_EQ(C.localLength(),2);
    EXPECT_DOUBLE_EQ( C(0,0), 3. );
    EXPECT_DOUBLE_EQ( C(0,1), 7. );
    EXPECT_DOUBLE_EQ( C(0,2), 0. );
    EXPECT_DOUBLE_EQ( C(0,3), 0. );

    EXPECT_DOUBLE_EQ( C(1,0), 3. );
    EXPECT_DOUBLE_EQ( C(1,1), 4. );
    EXPECT_DOUBLE_EQ( C(1,2), 0. );
    EXPECT_DOUBLE_EQ( C(1,3), 0. );
  }
}




// TEST_F(epetraSparseMatR7MultiVectorR9C4Fixture,
//        CrsMatrixTransposeMVecEpetraProduct){

//   using mymat_w_t = rompp::containers::Matrix<Epetra_CrsMatrix>;
//   using mymv_w_t = rompp::containers::MultiVector<Epetra_MultiVector>;  
//   ASSERT_TRUE(rompp::containers::meta::is_sparse_matrix_wrapper_epetra<mymat_w_t>::value);
//   ASSERT_TRUE(rompp::containers::meta::is_multi_vector_wrapper_epetra<mymv_w_t>::value);

//   // fill data
//   fillCrsMatrix();
//   fillMultiVector();
  
//   // we need to call fillcomplete considering rectangular matrix because
//   // the multivector has 9 rows, while matrix has 7 rows
//   // here we want to do: C = A^T * mv so: 
//   const auto & domainMap = *dataMapMV_; // = row map of sparse matrix
//   const auto & rangeMap = *dataMapMV_; // = row map of multivector
//   sm_->fillingIsCompleted(domainMap, rangeMap);

//   // mymv_w_t C( rangeMap, mv_->globalNumVectors() );
//   // sm_->data()->Multiply(true, *mv_->data(), *C.data()); 
//   // // compute prduct
//   auto C = rompp::containers::ops::product(*sm_, *mv_, true);
//   C.data()->Print(std::cout);
// }
