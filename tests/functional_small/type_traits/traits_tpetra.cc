#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include "pressio_type_traits.hpp"

TEST(tpetra, MVTraits)
{
  using namespace pressio;

  using ST = double;
  using LO = int;
  using GO = unsigned long;
  typedef Tpetra::MultiVector<>::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  using T = Tpetra::MultiVector<ST, LO, GO, NT>;
  using mytraits = pressio::Traits<T>;
  static_assert(pressio::is_multi_vector_tpetra<T>::value,"");

  ::testing::StaticAssertTypeEq<typename
  				mytraits::scalar_type, double>();

  ::testing::StaticAssertTypeEq<typename
  				mytraits::local_ordinal_type,int>();

  ::testing::StaticAssertTypeEq<typename
  				mytraits::global_ordinal_type,unsigned long>();

  ::testing::StaticAssertTypeEq<typename
  				mytraits::data_map_type, map_type>();

  ::testing::StaticAssertTypeEq<typename
  				mytraits::node_type, NT>();

  ::testing::StaticAssertTypeEq<typename
  				mytraits::dot_type, double>();

  ::testing::StaticAssertTypeEq<typename
  				mytraits::mag_type, double>();

  ::testing::StaticAssertTypeEq<typename
  				mytraits::communicator_type,
  				Teuchos::RCP<const Teuchos::Comm<int>>
  				>();

  ASSERT_TRUE(mytraits::rank == 2);
  ASSERT_TRUE(mytraits::is_shared_mem == false);

  ASSERT_TRUE(mytraits::multi_vector_identifier 
      == pressio::MultiVectorIdentifier::Tpetra);
}


TEST(tpetra, VectorTraits)
{
  using namespace pressio;

  using ST = double;
  using LO = int;
  using GO = unsigned long;
  typedef Tpetra::Vector<>::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  
  using T = Tpetra::Vector<ST, LO, GO, NT>;
  static_assert(::pressio::is_vector_tpetra<T>::value,"");
  using mytraits = pressio::Traits<T>;
 
  ::testing::StaticAssertTypeEq<typename
          mytraits::scalar_type, double>();

  ::testing::StaticAssertTypeEq<typename
          mytraits::local_ordinal_type,int>();

  ::testing::StaticAssertTypeEq<typename
          mytraits::global_ordinal_type,unsigned long>();

  ::testing::StaticAssertTypeEq<typename
          mytraits::data_map_type, map_type>();
  
  ::testing::StaticAssertTypeEq<typename
        mytraits::node_type, NT>();

  ::testing::StaticAssertTypeEq<typename
          mytraits::dot_type, double>();

  ::testing::StaticAssertTypeEq<typename
          mytraits::mag_type, double>();  

  ::testing::StaticAssertTypeEq<typename
          mytraits::communicator_type,
        Teuchos::RCP<const Teuchos::Comm<int>>
        >();
  
  ASSERT_TRUE(mytraits::rank == 1);
  ASSERT_TRUE(mytraits::is_shared_mem == false);
  ASSERT_TRUE(mytraits::is_dynamic == true);

  ASSERT_TRUE(mytraits::vector_identifier
        == pressio::VectorIdentifier::Tpetra);
}
