#include <gtest/gtest.h>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Version.h"
#include "Epetra_MpiComm.h"
#include "pressio_type_traits.hpp"

TEST(epetra, VectorTraits)
{
  using namespace pressio;

  using T = Epetra_Vector;
  static_assert(pressio::is_vector_epetra<T>::value,"");

  using vecTrait = pressio::traits<T>;
 
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::scalar_type, double>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::local_ordinal_type, int>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::global_ordinal_type, int>();
  
  ::testing::StaticAssertTypeEq<typename
  				vecTrait::data_map_type, Epetra_BlockMap>();

  ::testing::StaticAssertTypeEq<typename
  				vecTrait::communicator_type, Epetra_Comm>();
  
  ASSERT_TRUE(vecTrait::rank == 1);
  ASSERT_TRUE(vecTrait::is_shared_mem == false);
}

TEST(eped_epetra, MVTraits)
{
  using namespace pressio;

  using T = Epetra_MultiVector;
  static_assert(::pressio::is_multi_vector_epetra<T>::value, "");

  using vecTrait = pressio::traits<T>;
  ::testing::StaticAssertTypeEq<typename
       vecTrait::scalar_type, double>();

  ::testing::StaticAssertTypeEq<typename
       vecTrait::local_ordinal_type, int>();

  ::testing::StaticAssertTypeEq<typename
       vecTrait::global_ordinal_type, int>();

  ::testing::StaticAssertTypeEq<typename
       vecTrait::data_map_type,
       Epetra_BlockMap>();

  ::testing::StaticAssertTypeEq<typename
       vecTrait::communicator_type,
       Epetra_Comm>();
  
  ASSERT_TRUE(vecTrait::rank == 2);
  ASSERT_TRUE(vecTrait::is_shared_mem == false);
  ASSERT_TRUE(vecTrait::package_identifier == pressio::PackageIdentifier::Trilinos);
}