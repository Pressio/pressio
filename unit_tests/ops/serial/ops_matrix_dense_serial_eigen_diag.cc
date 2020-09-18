
#include <gtest/gtest.h>
#include "pressio_ops.hpp"
namespace{

  template <typename matrix_t, typename vec_t>
  void testAddition(matrix_t A, vec_t v, vec_t v2)
  {
    {
      auto diagvals = pressio::containers::diag(A);
      ::pressio::ops::do_update(diagvals,1.,v,2.);
      EXPECT_DOUBLE_EQ( diagvals(0), 3.2 );
      EXPECT_DOUBLE_EQ( diagvals(1), 10.2 );
      EXPECT_DOUBLE_EQ( diagvals(2), 17.2 );
      EXPECT_DOUBLE_EQ( diagvals(3), 24.  );
      ::pressio::ops::do_update(diagvals,v,2.);
      EXPECT_DOUBLE_EQ( diagvals(0), 2. );
      EXPECT_DOUBLE_EQ( diagvals(1), 4. );
      EXPECT_DOUBLE_EQ( diagvals(2), 6. );
      EXPECT_DOUBLE_EQ( diagvals(3), 8. );
      ::pressio::ops::do_update(diagvals,2.,v,2.,v2,2.);
      EXPECT_DOUBLE_EQ( diagvals(0), 28. );
      EXPECT_DOUBLE_EQ( diagvals(1), 36.  );
      ::pressio::ops::do_update(diagvals,v,2.,v2,2.);
      EXPECT_DOUBLE_EQ( diagvals(0), 24. );
      EXPECT_DOUBLE_EQ( diagvals(1), 28. );

      //cycle through rest of ops to check
      ::pressio::ops::do_update(diagvals,4.,v,2.,v2,2.,v2,2.);
      ::pressio::ops::do_update(diagvals,v,2.,v2,2.,v2,2.);
      ::pressio::ops::do_update(diagvals,4.,v,2.,v2,2.,v2,2.,v2,3.);
      ::pressio::ops::do_update(diagvals,v,2.,v2,2.,v2,2.,v2,3.);
    }
  }

  template <typename matrix_t, typename vec_t>
  void testDot(matrix_t A, vec_t v)
  {
    {
      auto diagvals = pressio::containers::diag(A);
      auto result = ::pressio::ops::dot(diagvals,v);
      EXPECT_NEAR( result, 111.2,1e-10);
      result = 0;
      ::pressio::ops::dot(diagvals,v,result);
      EXPECT_NEAR( result, 111.2,1e-10);
    }
  }


  template <typename matrix_t, typename vec_t>
  void testProduct(matrix_t A, vec_t result)
  {
    {
      auto diagvals = pressio::containers::diag(A);
      ::pressio::ops::product(::pressio::nontranspose(),
                              1.,A,diagvals,0.,result);
      EXPECT_NEAR( result(0), 111.44,1e-10);
      EXPECT_NEAR( result(1), 250.84,1e-10);
      EXPECT_NEAR( result(2), 390.24,1e-10);
      EXPECT_NEAR( result(3), 526.4 ,1e-10);
      ::pressio::ops::product(::pressio::transpose(),
                              1.,A,diagvals,0.,result);
      EXPECT_NEAR( result(0), 341.24,1e-10);
      EXPECT_NEAR( result(1), 376.84,1e-10);
      EXPECT_NEAR( result(2), 412.44,1e-10);
      EXPECT_NEAR( result(3), 444.8,1e-10);
    }
  }


};


TEST(containers_matrix_serial_eigen_diag_ops, diag)
{
  // col-major matrix (which is default in Eigen)
  using eigmat_t = Eigen::Matrix<double,-1,-1>;
  using eigmat_t_float = Eigen::Matrix<float,-1,-1>;
  using eigvec_t = Eigen::Matrix<double, -1, 1>;

  using myM_t = pressio::containers::DenseMatrix<eigmat_t>;
  using myM_t_float = pressio::containers::DenseMatrix<eigmat_t_float>;
  using myV_t = pressio::containers::Vector<eigvec_t>;


  myM_t A(4,4);
  myM_t_float Af(4,4);

  myV_t v1(4);
  myV_t v2(4);

  A(0,0) = 1.2;  A(0,1) = 2.;  A(0,2) = 3.;  A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.2;  A(1,2) = 7.;  A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.2; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.; A(3,3) = 16.;
  
  v1(0) = 1.;
  v1(1) = 2.;
  v1(2) = 3.;
  v1(3) = 4.;

  v2(0) = 11.;
  v2(1) = 12.;
  v2(2) = 13.;
  v2(3) = 14.;

  testAddition(A,v1,v2);
  testDot(A,v1);
  testProduct(A,v1);

}

TEST(containers_matrix_serial_eigen_diag_ops, diagRowMajor)
{
  // col-major matrix (which is default in Eigen)
  using eigmat_t = Eigen::Matrix<double,-1,-1,Eigen::RowMajor>;
  using eigvec_t = Eigen::Matrix<double, -1, 1>;

  using myM_t = pressio::containers::DenseMatrix<eigmat_t>;
  using myV_t = pressio::containers::Vector<eigvec_t>;


  myM_t A(4,4);

  myV_t v1(4);
  myV_t v2(4);

  A(0,0) = 1.2;  A(0,1) = 2.;  A(0,2) = 3.;  A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.2;  A(1,2) = 7.;  A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.2; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.; A(3,3) = 16.;
  
  v1(0) = 1.;
  v1(1) = 2.;
  v1(2) = 3.;
  v1(3) = 4.;

  v2(0) = 11.;
  v2(1) = 12.;
  v2(2) = 13.;
  v2(3) = 14.;

  testAddition(A,v1,v2);
  testDot(A,v1);
  testProduct(A,v1);
}
