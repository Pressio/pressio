#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_SerialDenseVector.hpp"
#include <Teuchos_SerialDenseMatrix.hpp>
#include "pressio/type_traits.hpp"

TEST(type_traits, isTeuchosRCP)
{
  using namespace pressio;

  class foo{
    int a_ = 0;
    public:
      foo(int a) : a_(a) {};
  };

  using foo_t1 = foo;
  using foo_t2 = foo *;
  using foo_t3 = std::shared_ptr<foo>;
  using foo_t4 = Teuchos::RCP<foo>;
  using foo_t5 = Teuchos::RCP<const foo>;

  EXPECT_EQ( pressio::is_teuchos_rcp<foo_t1>::value, false);
  EXPECT_EQ( pressio::is_teuchos_rcp<foo_t2>::value, false);
  EXPECT_EQ( pressio::is_teuchos_rcp<foo_t3>::value, false);
  EXPECT_EQ( pressio::is_teuchos_rcp<foo_t4>::value, true);
  EXPECT_EQ( pressio::is_teuchos_rcp<foo_t5>::value, true);
}

TEST(type_traits, TeuchosVector)
{
  using namespace pressio;

  using T = Teuchos::SerialDenseVector<int, double>;
  using mytraits = pressio::Traits<T>;
  static_assert(pressio::is_dense_vector_teuchos<T>::value,"");

  ::testing::StaticAssertTypeEq<typename mytraits::scalar_type, double>();

  ::testing::StaticAssertTypeEq<typename mytraits::ordinal_type,int>();

  ASSERT_TRUE(mytraits::rank == 1);
  ASSERT_TRUE(mytraits::is_shared_mem == true);

  ASSERT_TRUE(mytraits::vector_identifier 
      == pressio::VectorIdentifier::TeuchosSerialDense);
}

TEST(type_traits, TeuchosMatrix)
{
  using namespace pressio;

  using T = Teuchos::SerialDenseMatrix<long long, double>;
  using mytraits = pressio::Traits<T>;
  static_assert(pressio::is_dense_matrix_teuchos<T>::value,"");

  ::testing::StaticAssertTypeEq<typename mytraits::scalar_type, double>();
  ::testing::StaticAssertTypeEq<typename mytraits::ordinal_type, long long>();

  ASSERT_TRUE(mytraits::rank == 2);
  ASSERT_TRUE(mytraits::is_shared_mem == true);

  ASSERT_TRUE(mytraits::matrix_identifier 
      == pressio::MatrixIdentifier::DenseTeuchosSerial);
}
