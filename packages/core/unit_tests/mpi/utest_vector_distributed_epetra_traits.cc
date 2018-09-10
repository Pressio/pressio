#include <gtest/gtest.h>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#include "Epetra_MpiComm.h"
#include "CORE_VECTOR"

TEST(core_vector_distributed_epetra, EpetraVectorTraits)
{
  using natV_t = Epetra_Vector;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(natV_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(natV_t);
  STATIC_ASSERT_IS_VECTOR_EPETRA(natV_t);

  using myvec_t = core::Vector<natV_t>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(myvec_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(myvec_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(myvec_t);

  using vecTrait = core::details::traits<myvec_t>;
 
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::scalar_t,
  				core::defaultTypes::epetra_scalar_t>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::local_ordinal_t,
  				core::defaultTypes::epetra_lo_t>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::global_ordinal_t,
  				core::defaultTypes::epetra_go_t1>();
  
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::wrapped_t, natV_t>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::derived_t, myvec_t>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::data_map_t, Epetra_BlockMap>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::communicator_t, Epetra_Comm>();
  
  ASSERT_TRUE(vecTrait::is_vector == 1);
  ASSERT_TRUE(vecTrait::is_shared_mem == 0);

}//end TEST
