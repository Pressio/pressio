
#include <gtest/gtest.h>
// #include "matrix/concrete/core_matrix_sparse_serial_eigen.hpp"
// #include "matrix/meta/core_matrix_meta.hpp"
#include "CORE_MATRIX"


template<typename scalar, int storage, typename index>
struct typesEig{
  using sc_t = scalar;
  using mat_t = Eigen::SparseMatrix<sc_t,storage,index>;
  using storeindex_t = index;
  static constexpr int storemode = storage;
};


template<typename T>
struct core_matrix_sparse_serial_eigen_classTest
  : public ::testing::Test{
public:
  using native_t = typename T::mat_t;
  using mymat_t = core::matrix<native_t>;
  using matTrait = core::details::traits<mymat_t>;

  mymat_t createMatrix1(){
    mymat_t m4(6,6);
    using veci = std::vector<int>;
    using vecd = std::vector<double>;
    std::map<int,std::pair<vecd,veci>> DD;
    DD[0] = std::make_pair<vecd,veci>({2.,4.},   {1,2});
    DD[2] = std::make_pair<vecd,veci>({3.,5.,6.},{3,4,5});
    DD[5] = std::make_pair<vecd,veci>({1.,22.},  {2,4});
    for (auto const & it : DD){
      int ind = it.first;
      vecd data = std::get<0>(it.second);
      veci cols = std::get<1>(it.second);
      m4.insertValues(ind,(int)data.size(),data.data(),cols.data());
    }
    return m4;
  };
  mymat_t createMatrix2(){
    mymat_t m4(6,6);
    using veci = std::vector<int>;
    using vecd = std::vector<double>;
    std::map<int,std::pair<vecd,veci>> DD;
    DD[0] = std::make_pair<vecd,veci>({5.,7.},   {1,2});
    DD[1] = std::make_pair<vecd,veci>({2.,15.,6.},{2,3,4});
    DD[5] = std::make_pair<vecd,veci>({10.,22.},  {2,4});
    for (auto const & it : DD){
      int ind = it.first;
      vecd data = std::get<0>(it.second);
      veci cols = std::get<1>(it.second);
      m4.insertValues(ind,(int)data.size(),data.data(),cols.data());
    }
    return m4;
  };

  void checkAdditionOperator(){
    mymat_t m1 = createMatrix1();
    mymat_t m2 = createMatrix2();
    mymat_t m3 = m1 + m2;

    // the inputs have the same dimensions
    EXPECT_EQ( m1.nonZerosCount(), m2.nonZerosCount() );
    // the result does not necessarily because it is sum
    EXPECT_NE( m3.nonZerosCount(), m1.nonZerosCount() );

    if (T::storemode==Eigen::RowMajor){
      EXPECT_DOUBLE_EQ( m3(0,1), 7. );
      EXPECT_DOUBLE_EQ( m3(0,2), 11. );
      EXPECT_DOUBLE_EQ( m3(1,2), 2. );
      EXPECT_DOUBLE_EQ( m3(1,3), 15. );
      EXPECT_DOUBLE_EQ( m3(1,4), 6. );
      EXPECT_DOUBLE_EQ( m3(2,3), 3. );
      EXPECT_DOUBLE_EQ( m3(2,4), 5. );
      EXPECT_DOUBLE_EQ( m3(2,5), 6. );
      EXPECT_DOUBLE_EQ( m3(5,2), 11. );
      EXPECT_DOUBLE_EQ( m3(5,3), 0. );
      EXPECT_DOUBLE_EQ( m3(5,4), 44. );
    }
    if (T::storemode == Eigen::ColMajor){
      EXPECT_DOUBLE_EQ( m3(1,0), 7. );
      EXPECT_DOUBLE_EQ( m3(2,0), 11. );
      EXPECT_DOUBLE_EQ( m3(2,1), 2. );
      EXPECT_DOUBLE_EQ( m3(3,1), 15. );
      EXPECT_DOUBLE_EQ( m3(4,1), 6. );
      EXPECT_DOUBLE_EQ( m3(3,2), 3. );
      EXPECT_DOUBLE_EQ( m3(4,2), 5. );
      EXPECT_DOUBLE_EQ( m3(5,2), 6. );
      EXPECT_DOUBLE_EQ( m3(2,5), 11. );
      EXPECT_DOUBLE_EQ( m3(3,5), 0. );
      EXPECT_DOUBLE_EQ( m3(4,5), 44. );
    }
  }
  //---------------------------------------------------

  void checkCompoundAddOperator(){
    mymat_t m1 = createMatrix1();
    mymat_t m2 = createMatrix2();
    // the inputs have the same dimensions
    EXPECT_EQ( m1.nonZerosCount(), m2.nonZerosCount() );
    m1 += m2;
    // the result does not necessarily because it is sum
    EXPECT_NE( m1.nonZerosCount(), m2.nonZerosCount() );

    if (T::storemode==Eigen::RowMajor){
      EXPECT_DOUBLE_EQ( m1(0,1), 7. );
      EXPECT_DOUBLE_EQ( m1(0,2), 11. );
      EXPECT_DOUBLE_EQ( m1(1,2), 2. );
      EXPECT_DOUBLE_EQ( m1(1,3), 15. );
      EXPECT_DOUBLE_EQ( m1(1,4), 6. );
      EXPECT_DOUBLE_EQ( m1(2,3), 3. );
      EXPECT_DOUBLE_EQ( m1(2,4), 5. );
      EXPECT_DOUBLE_EQ( m1(2,5), 6. );
      EXPECT_DOUBLE_EQ( m1(5,2), 11. );
      EXPECT_DOUBLE_EQ( m1(5,3), 0. );
      EXPECT_DOUBLE_EQ( m1(5,4), 44. );
    }
    if (T::storemode == Eigen::ColMajor){
      EXPECT_DOUBLE_EQ( m1(1,0), 7. );
      EXPECT_DOUBLE_EQ( m1(2,0), 11. );
      EXPECT_DOUBLE_EQ( m1(2,1), 2. );
      EXPECT_DOUBLE_EQ( m1(3,1), 15. );
      EXPECT_DOUBLE_EQ( m1(4,1), 6. );
      EXPECT_DOUBLE_EQ( m1(3,2), 3. );
      EXPECT_DOUBLE_EQ( m1(4,2), 5. );
      EXPECT_DOUBLE_EQ( m1(5,2), 6. );
      EXPECT_DOUBLE_EQ( m1(2,5), 11. );
      EXPECT_DOUBLE_EQ( m1(3,5), 0. );
      EXPECT_DOUBLE_EQ( m1(4,5), 44. );
    }
  }
  //---------------------------------------------------

  
  void checkSubtractionOperator(){
    mymat_t m1 = createMatrix1();
    mymat_t m2 = createMatrix2();
    mymat_t m3 = m1 - m2;
    // the inputs have the same dimensions
    EXPECT_EQ( m1.nonZerosCount(), m2.nonZerosCount() );
    // the result does not necessarily because it is sum
    EXPECT_NE( m3.nonZerosCount(), m1.nonZerosCount() );

    if (T::storemode==Eigen::RowMajor){
      EXPECT_DOUBLE_EQ( m3(0,1), -3. );
      EXPECT_DOUBLE_EQ( m3(0,2), -3. );
      EXPECT_DOUBLE_EQ( m3(1,2), -2. );
      EXPECT_DOUBLE_EQ( m3(1,3), -15. );
      EXPECT_DOUBLE_EQ( m3(1,4), -6. );
      EXPECT_DOUBLE_EQ( m3(2,3), 3. );
      EXPECT_DOUBLE_EQ( m3(2,4), 5. );
      EXPECT_DOUBLE_EQ( m3(2,5), 6. );
      EXPECT_DOUBLE_EQ( m3(5,2), -9 );
      EXPECT_DOUBLE_EQ( m3(5,3), 0. );
      EXPECT_DOUBLE_EQ( m3(5,4), 0. );
    }
    if (T::storemode == Eigen::ColMajor){
      EXPECT_DOUBLE_EQ( m3(1,0), -3. );
      EXPECT_DOUBLE_EQ( m3(2,0), -3. );
      EXPECT_DOUBLE_EQ( m3(2,1), -2. );
      EXPECT_DOUBLE_EQ( m3(3,1), -15. );
      EXPECT_DOUBLE_EQ( m3(4,1), -6. );
      EXPECT_DOUBLE_EQ( m3(3,2), 3. );
      EXPECT_DOUBLE_EQ( m3(4,2), 5. );
      EXPECT_DOUBLE_EQ( m3(5,2), 6. );
      EXPECT_DOUBLE_EQ( m3(2,5), -9. );
      EXPECT_DOUBLE_EQ( m3(3,5), 0. );
      EXPECT_DOUBLE_EQ( m3(4,5), 0. );
    }
  }
  //---------------------------------------------------

  void checkCompoundSubOperator(){
    mymat_t m1 = createMatrix1();
    mymat_t m2 = createMatrix2();
    // the inputs have the same dimensions
    EXPECT_EQ( m1.nonZerosCount(), m2.nonZerosCount() );
    m1 -= m2;
    // the result does not necessarily because it is sum
    EXPECT_NE( m1.nonZerosCount(), m2.nonZerosCount() );

    if (T::storemode==Eigen::RowMajor){
      EXPECT_DOUBLE_EQ( m1(0,1), -3. );
      EXPECT_DOUBLE_EQ( m1(0,2), -3. );
      EXPECT_DOUBLE_EQ( m1(1,2), -2. );
      EXPECT_DOUBLE_EQ( m1(1,3), -15. );
      EXPECT_DOUBLE_EQ( m1(1,4), -6. );
      EXPECT_DOUBLE_EQ( m1(2,3), 3. );
      EXPECT_DOUBLE_EQ( m1(2,4), 5. );
      EXPECT_DOUBLE_EQ( m1(2,5), 6. );
      EXPECT_DOUBLE_EQ( m1(5,2), -9 );
      EXPECT_DOUBLE_EQ( m1(5,3), 0. );
      EXPECT_DOUBLE_EQ( m1(5,4), 0. );
    }
    if (T::storemode == Eigen::ColMajor){
      EXPECT_DOUBLE_EQ( m1(1,0), -3. );
      EXPECT_DOUBLE_EQ( m1(2,0), -3. );
      EXPECT_DOUBLE_EQ( m1(2,1), -2. );
      EXPECT_DOUBLE_EQ( m1(3,1), -15. );
      EXPECT_DOUBLE_EQ( m1(4,1), -6. );
      EXPECT_DOUBLE_EQ( m1(3,2), 3. );
      EXPECT_DOUBLE_EQ( m1(4,2), 5. );
      EXPECT_DOUBLE_EQ( m1(5,2), 6. );
      EXPECT_DOUBLE_EQ( m1(2,5), -9. );
      EXPECT_DOUBLE_EQ( m1(3,5), 0. );
      EXPECT_DOUBLE_EQ( m1(4,5), 0. );
    }
  }
  //---------------------------------------------------
  

  void checkStarOperator(){
    mymat_t m1 = createMatrix1();
    mymat_t m2 = createMatrix2();
    mymat_t m3 = m1 * m2;
    // the inputs have the same dimensions
    EXPECT_EQ( m1.nonZerosCount(), m2.nonZerosCount() );
    // the result does not necessarily because it is sum
    EXPECT_NE( m3.nonZerosCount(), m1.nonZerosCount() );

    if (T::storemode==Eigen::RowMajor){
      EXPECT_DOUBLE_EQ( m3(0,1), 0. );
      EXPECT_DOUBLE_EQ( m3(0,2), 4. );
      EXPECT_DOUBLE_EQ( m3(0,3), 30. );
    }
    if (T::storemode == Eigen::ColMajor){
      EXPECT_DOUBLE_EQ( m3(4,0), 35. );
    }
  }
  //---------------------------------------------------
  

  //---------------------------------------------------
  void checkConstructor(){
    STATIC_ASSERT_IS_MATRIX_SPARSE_SERIAL_EIGEN(native_t);
    EXPECT_TRUE(matTrait::isEigen == 1);
    //EXPECT_TRUE(matTrait::isRowMajor == storage_type);
    // mymat_t m1;
    // EXPECT_EQ( m1.rows(), 0 );
    // EXPECT_EQ( m1.cols(), 0 );
    mymat_t m2(5, 8);
    EXPECT_EQ( m2.rows(), 5 );
    EXPECT_EQ( m2.cols(), 8 );
    native_t eigMat(4,5);
    mymat_t m3(eigMat);
    EXPECT_EQ( m3.rows(), 4 );
    EXPECT_EQ( m3.cols(), 5 );
  };//end constr
  //---------------------------------------------------

  void checkFilling(){
    // create matrix 
    mymat_t m4(12,12);
    // fill
    using veci = std::vector<int>;
    using vecd = std::vector<double>;
    std::map<int,std::pair<vecd,veci>> DD;
    DD[0] = std::make_pair<vecd,veci>({2.3,4.1}, {1,2});
    DD[3] = std::make_pair<vecd,veci>({2.3,4.1,6.6}, {5,4,6});
    DD[7] = std::make_pair<vecd,veci>({1.,3.}, {7,9});
    for (auto const & it : DD){
      int ind = it.first;
      vecd data = std::get<0>(it.second);
      veci cols = std::get<1>(it.second);
      m4.insertValues(ind,(int)data.size(),data.data(),cols.data());
    }
    // check
    EXPECT_EQ( m4.nonZerosCount(), 7 );
    if (T::storemode==Eigen::RowMajor){
      EXPECT_DOUBLE_EQ( m4(0,0), 0.0 );
      EXPECT_DOUBLE_EQ( m4(0,1), 2.3 );
      EXPECT_DOUBLE_EQ( m4(0,2), 4.1 );
      EXPECT_DOUBLE_EQ( m4(3,5), 2.3 );
      EXPECT_DOUBLE_EQ( m4(3,4), 4.1 );
      EXPECT_DOUBLE_EQ( m4(3,6), 6.6 );
      EXPECT_DOUBLE_EQ( m4(7,7), 1. );
      EXPECT_DOUBLE_EQ( m4(7,9), 3. );
      EXPECT_DOUBLE_EQ( m4(7,10), 0.0 );
    }
    if (T::storemode == Eigen::ColMajor){
      EXPECT_DOUBLE_EQ( m4(0,0), 0.0 );
      EXPECT_DOUBLE_EQ( m4(1,0), 2.3 );
      EXPECT_DOUBLE_EQ( m4(2,0), 4.1 );
      EXPECT_DOUBLE_EQ( m4(5,3), 2.3 );
      EXPECT_DOUBLE_EQ( m4(4,3), 4.1 );
      EXPECT_DOUBLE_EQ( m4(6,3), 6.6 );
      EXPECT_DOUBLE_EQ( m4(7,7), 1. );
      EXPECT_DOUBLE_EQ( m4(9,7), 3. );
      EXPECT_DOUBLE_EQ( m4(10,7), 0.0 );
    }
  };//end 
  //---------------------------------------------------

  void checkQuery(){
    // mymat_t m1;
    // ::testing::StaticAssertTypeEq<decltype(m1.data()),
				//   native_t * >();
    const mymat_t m2(45,64);
    ::testing::StaticAssertTypeEq< decltype(m2.data()),
				   const native_t * >();
  };
  //---------------------------------------------------
};


typedef ::testing::Types< typesEig<double,Eigen::RowMajor,int>,
			  typesEig<double,Eigen::ColMajor,int>
			  > MyTypes;
TYPED_TEST_CASE(core_matrix_sparse_serial_eigen_classTest, MyTypes);

TYPED_TEST(core_matrix_sparse_serial_eigen_classTest, constructor){
  this->checkConstructor();
}

TYPED_TEST(core_matrix_sparse_serial_eigen_classTest, filling){
  this->checkFilling();
}

TYPED_TEST(core_matrix_sparse_serial_eigen_classTest, queryWrappedData){
  this->checkQuery();
}

TYPED_TEST(core_matrix_sparse_serial_eigen_classTest, additionOperator){
  this->checkAdditionOperator();
}

TYPED_TEST(core_matrix_sparse_serial_eigen_classTest, subtractionOperator){
  this->checkSubtractionOperator();
}

TYPED_TEST(core_matrix_sparse_serial_eigen_classTest, starOperator){
  this->checkStarOperator();
}

TYPED_TEST(core_matrix_sparse_serial_eigen_classTest, compoundAddOperator){
  this->checkCompoundAddOperator();
}

TYPED_TEST(core_matrix_sparse_serial_eigen_classTest, compoundSubOperator){
  this->checkCompoundSubOperator();
}


