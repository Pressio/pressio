#include <gtest/gtest.h>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#include "Epetra_MpiComm.h"
#include "pressio_containers.hpp"

TEST(containers_vector_distributed_epetra, EpetraVectorTraits)
{
  using namespace pressio;

  using natV_t = Epetra_Vector;
  using myvec_t = containers::Vector<natV_t>;

  static_assert(::pressio::containers::predicates::is_vector_epetra<natV_t>::value,"");
  static_assert(::pressio::containers::predicates::is_vector_wrapper_epetra<myvec_t>::value,"");

  using vecTrait = containers::details::traits<myvec_t>;
 
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::scalar_t, double>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::local_ordinal_t, int>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::global_ordinal_t, int>();
  
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::wrapped_t, natV_t>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::derived_t, myvec_t>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::data_map_t, Epetra_BlockMap>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::communicator_t, Epetra_Comm>();
  
  ASSERT_TRUE(vecTrait::is_vector == true);
  ASSERT_TRUE(vecTrait::is_shared_mem == false);

}//end TEST
