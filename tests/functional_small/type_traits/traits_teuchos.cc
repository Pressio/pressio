#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "traits_shared.hpp"

namespace pressio { namespace traits { namespace test {

template <typename T, int rank>
void test_teuchos_container()
{
  // traits and shared predicates
  test_container_traits<
    T,
    rank,
    typename T::scalarType
  >();

  // negative checks (cross-package)
  test_is_not_eigen_container<T>();
  test_is_not_kokkos_container<T>();
  test_is_not_tpetra_container<T>();
  test_is_not_tpetra_block_container<T>();
  test_is_not_epetra_container<T>();
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

  // traits and shared predicates
  test_teuchos_container<T, 1>();

  // vector predicates
  static_assert(pressio::is_dense_vector_teuchos<T>::value,"");

  // negative checks (within Teuchos)
  static_assert(pressio::is_dense_matrix_teuchos<T>::value == false, "");
}

TEST(type_traits, TeuchosMatrix)
{
  using T = Teuchos::SerialDenseMatrix<long long, float>;

  // traits and shared predicates
  test_teuchos_container<T, 2>();

  // matrix predicates
  static_assert(pressio::is_dense_matrix_teuchos<T>::value,"");

  // negative checks (within Teuchos)
  static_assert(pressio::is_dense_vector_teuchos<T>::value == false, "");
}

}}} // pressio::traits::test
