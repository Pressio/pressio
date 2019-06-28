#include <gtest/gtest.h>
#include "ALGEBRA_VECTOR"
#include "ALGEBRA_MULTI_VECTOR"
#include "ALGEBRA_MATRIX"

TEST(algebra_matrxi_dense_distributed_epetra, EpetraMultiVectorActingDenseMatrixTraits)
{
  using natV_t = Epetra_MultiVector;

  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(natV_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(natV_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(natV_t);

  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SHAREDMEM_EIGEN(natV_t);
  STATIC_ASSERT_IS_NOT_MATRIX_SPARSE_SHAREDMEM_EIGEN(natV_t);
  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SHAREDMEM_STDLIB(natV_t);
  STATIC_ASSERT_IS_NOT_MATRIX_SPARSE_DISTRIBUTED_EPETRA(natV_t);

  // the native multivector type itself can be seen as a multiver or dense matrix
  STATIC_ASSERT_IS_MULTIVECTOR_EPETRA(natV_t);
  STATIC_ASSERT_IS_MATRIX_DENSE_DISTRIBUTED_EPETRA(natV_t);
  // if passed to a Matrix wrapper then it becomes a distr dense matrix
  // if passed to a multivector, then it becomes a distr multivector
  using myDDM_t = algebra::Matrix<natV_t>;
  STATIC_ASSERT_IS_ALGEBRA_MATRIX_WRAPPER(myDDM_t);
  ASSERT_FALSE( algebra::meta::is_algebra_multi_vector_wrapper<myDDM_t>::value );
  //  STATIC_ASSERT_IS_NOT_ALGEBRA_MULTI_VECTOR_WRAPPER(myDDM_t);

  using myDMV_t = algebra::MultiVector<natV_t>;
  STATIC_ASSERT_IS_ALGEBRA_MULTI_VECTOR_WRAPPER(myDMV_t);
  // STATIC_ASSERT_IS_NOT_ALGEBRA_MATRIX_WRAPPER(myDMV_t);

  
  using mytraits = algebra::details::traits<myDDM_t>;
  ::testing::StaticAssertTypeEq<typename
  				mytraits::scalar_t,
  				algebra::default_types::epetra_scalar_t>();
  
  ASSERT_TRUE(mytraits::is_vector == 0);
  ASSERT_TRUE(mytraits::isEpetra == 1);
  ASSERT_TRUE(mytraits::actingAsMultiVector == 0);
  ASSERT_TRUE(mytraits::actingAsDenseMatrix == 1);
  ASSERT_TRUE(mytraits::isEigen == 0);
  ASSERT_TRUE(mytraits::is_shared_mem == 0);

}//end TEST


