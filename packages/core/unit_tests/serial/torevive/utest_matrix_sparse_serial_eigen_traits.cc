
#include <gtest/gtest.h>
#include "CORE_MATRIX"


template<typename scalar, int storage, typename index>
struct typesEig{
  using sc_t = scalar;
  using mat_t = Eigen::SparseMatrix<sc_t,storage,index>;
  using storeindex_t = index;
  static constexpr int storemode = storage;
};

template<typename T>
struct core_matrix_sparse_serial_eigen_traitsTest
  : public ::testing::Test{
public:
  using native_t = typename T::mat_t;
  STATIC_ASSERT_IS_NOT_MATRIX_DENSE_SHAREDMEM_EIGEN(native_t);
  STATIC_ASSERT_IS_MATRIX_SPARSE_SHAREDMEM_EIGEN(native_t);
  using my_t = core::Matrix<native_t>;
  using myTrait = core::details::traits<my_t>;

  //need this void method otherwise the static assertion for type dont work
  void check(){
    ::testing::StaticAssertTypeEq<typename myTrait::scalar_t,
				  typename T::sc_t>();
    ::testing::StaticAssertTypeEq<typename myTrait::ordinal_t,
				  typename T::storeindex_t>();
    ::testing::StaticAssertTypeEq<typename
				  myTrait::wrapped_t,native_t>();
    ::testing::StaticAssertTypeEq<typename
				  myTrait::derived_t,my_t>();

    if (T::storemode == Eigen::RowMajor)
      EXPECT_TRUE(myTrait::is_row_major == Eigen::StorageOptions::RowMajor);
    else
       EXPECT_TRUE(myTrait::is_row_major == Eigen::ColMajor);
    
    ASSERT_TRUE(myTrait::is_col_major == !myTrait::is_row_major);
    ASSERT_TRUE(myTrait::wrapped_matrix_identifier 
      == core::details::WrappedMatrixIdentifier::SparseEigen);
    ASSERT_TRUE(myTrait::is_dense == 0);
    ASSERT_TRUE(myTrait::is_sparse == 1);
    ASSERT_TRUE(myTrait::is_vector == 0);
    ASSERT_TRUE(myTrait::is_shared_mem == 1);
    ASSERT_TRUE(myTrait::is_static == 0);
  }
};

typedef ::testing::Types< typesEig<int,Eigen::RowMajor,int>,
			  typesEig<int,Eigen::ColMajor,short>,
			  typesEig<int,Eigen::RowMajor,long long>,
			  //
			  typesEig<double,Eigen::RowMajor,int>,
			  typesEig<double,Eigen::ColMajor,short>,
			  typesEig<double,Eigen::RowMajor,long long>,
			  //
			  typesEig<float,Eigen::RowMajor,int>,
			  typesEig<float,Eigen::ColMajor,short>,
			  typesEig<float,Eigen::RowMajor,long long>,
			  //
			  typesEig<unsigned int,Eigen::RowMajor,int>,
			  typesEig<unsigned int,Eigen::ColMajor,short>,
			  typesEig<unsigned int,Eigen::RowMajor,long long>,
			  //
			  typesEig<std::complex<double>,Eigen::RowMajor,int>,
			  typesEig<std::complex<double>,Eigen::ColMajor,short>,
			  typesEig<std::complex<double>,Eigen::RowMajor,long long>
			  > MyTypes;
TYPED_TEST_CASE(core_matrix_sparse_serial_eigen_traitsTest, MyTypes);

TYPED_TEST(core_matrix_sparse_serial_eigen_traitsTest, traits)
{
  //this runs all types, no need to put anything
  this->check();

  using type1 = core::Matrix<Eigen::SparseMatrix<double,Eigen::RowMajor,int>>;  
  using type2 = core::Matrix<Eigen::SparseMatrix<double,Eigen::ColMajor,int>>;
  static_assert( core::meta::sparse_sharedmem_eigen_same_storage<type1,type2>::value == false, "ohh" );
}
