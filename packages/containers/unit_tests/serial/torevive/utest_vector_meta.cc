#include <gtest/gtest.h>
#include "CONTAINERS_VECTOR"


template<typename scalar, int rows, int cols>
struct typesEig{
  using sc_t = scalar;
  constexpr static int nr = rows;
  constexpr static int nc = cols;
};


template <typename T>
class VectorEigenMetaTest : public ::testing::Test {
	using namespace pressio;

	public:
  using native_t = Eigen::Matrix<typename T::sc_t, T::nr, T::nc>;

  STATIC_ASSERT_IS_VECTOR_EIGEN(native_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(native_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(native_t); 
};

typedef ::testing::Types< typesEig<int,-1,1>,
			  typesEig<int,5,1>,
			  typesEig<int,1,-1>,
			  typesEig<int,1,5>,
			  //
			  typesEig<double,-1,1>,
			  typesEig<double,5,1>,
			  typesEig<double,1,-1>,
			  typesEig<double,1,5>,
			  //
			  typesEig<float,-1,1>,
			  typesEig<float,5,1>,
			  typesEig<float,1,-1>,
			  typesEig<float,1,5>,
			  //
			  typesEig<std::complex<double>,-1,1>,
			  typesEig<std::complex<double>,5,1>,
			  typesEig<std::complex<double>,1,-1>,
			  typesEig<std::complex<double>,1,5>,
			  //
			  typesEig<unsigned int,-1,1>,
			  typesEig<unsigned int,5,1>,
			  typesEig<unsigned int,1,-1>,
			  typesEig<unsigned int,1,5>
			  > MyTypes;
TYPED_TEST_CASE(VectorEigenMetaTest, MyTypes);

TYPED_TEST(VectorEigenMetaTest, meta)
{
  //this runs all types, no need to put anything
}

TEST(containers_vector_meta, vectorMetasStdlib)
{
	using namespace pressio;

  using native_t = std::vector<double>;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(native_t);
  STATIC_ASSERT_IS_VECTOR_STDLIB(native_t);
  STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(native_t);
}

TEST(containers_vector_meta, vectorMetasEpetra)
{
	using namespace pressio;

  using native_t = Epetra_Vector;
  STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(native_t);
  STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(native_t);
  STATIC_ASSERT_IS_VECTOR_EPETRA(native_t);
}
