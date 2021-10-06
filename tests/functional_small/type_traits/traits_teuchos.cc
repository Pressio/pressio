#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "traits_shared.hpp"

template <typename T, int rank>
void test_teuchos_container()
{
  test_container_traits<
    T,
    pressio::PackageIdentifier::Trilinos,
    rank,
    true, /* shared mem */
    true, /* dynamic */
    typename T::scalarType,
    typename T::ordinalType
  >();
}

//*******************************
// Teuchos vector
//*******************************

TEST(type_traits, isTeuchosRCP)
{
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
  using T = Teuchos::SerialDenseVector<int, double>;
  static_assert(pressio::is_dense_vector_teuchos<T>::value,"");
  ASSERT_TRUE(pressio::Traits<T>::vector_identifier
      == pressio::VectorIdentifier::TeuchosSerialDense);
  test_teuchos_container<T, 1>();
}

TEST(type_traits, TeuchosMatrix)
{
  using T = Teuchos::SerialDenseMatrix<long long, float>;
  static_assert(pressio::is_dense_matrix_teuchos<T>::value,"");
  ASSERT_TRUE(pressio::Traits<T>::matrix_identifier
      == pressio::MatrixIdentifier::DenseTeuchosSerial);
  test_teuchos_container<T, 2>();
}
