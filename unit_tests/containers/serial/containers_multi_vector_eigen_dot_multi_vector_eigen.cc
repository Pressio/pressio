
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"

namespace{
  using eigdmat_t = Eigen::MatrixXd;
  using myMV_t = pressio::containers::MultiVector<eigdmat_t>;

  struct Data{
    myMV_t A;
    myMV_t B;

    Data() : A(6,3), B(6,4){
      A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
      A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
      A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
      A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
      A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
      A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;

      B(0,0) = 1.; B(0,1) = 2.; B(0,2) = 3.; B(0,3) = 3.;
      B(1,0) = 3.; B(1,1) = 2.; B(1,2) = 1.; B(1,3) = 3.;
      B(2,0) = 0.; B(2,1) = 0.; B(2,2) = 1.; B(2,3) = 3.;
      B(3,0) = 0.; B(3,1) = 1.; B(3,2) = 0.; B(3,3) = 3.;
      B(4,0) = 1.; B(4,1) = 0.; B(4,2) = 0.; B(4,3) = 3.;
      B(5,0) = 0.; B(5,1) = 1.; B(5,2) = 1.; B(5,3) = 3.;
    }
  };

  template <typename matrix_t>
  void doDot(const myMV_t & A, const myMV_t & B){
    auto C = pressio::containers::ops::dot<myMV_t, matrix_t>(A,B);
    ASSERT_EQ(C.rows(), 3);
    ASSERT_EQ(C.cols(), 4);

    // first row
    EXPECT_DOUBLE_EQ( C(0,0), 11.);
    EXPECT_DOUBLE_EQ( C(0,1), 8.);
    EXPECT_DOUBLE_EQ( C(0,2), 6.);
    EXPECT_DOUBLE_EQ( C(0,3), 15.);

    // last row
    EXPECT_DOUBLE_EQ( C(2,0), 6.);
    EXPECT_DOUBLE_EQ( C(2,1), 9.);
    EXPECT_DOUBLE_EQ( C(2,2), 12.);
    EXPECT_DOUBLE_EQ( C(2,3), 18.);
  }

  template <typename matrix_t>
  void doDot(const myMV_t & A, const myMV_t & B, matrix_t & C)
  {
    // subspan a block of size (3 x 4) so that we can store A^T B
    auto block = C.subspan(std::make_pair(10, 13), std::make_pair(15,19));
    pressio::containers::ops::dot(A,B,block);
    ASSERT_EQ(block.rows(), 3);
    ASSERT_EQ(block.cols(), 4);
    EXPECT_DOUBLE_EQ( block(0,0), 11.);
    EXPECT_DOUBLE_EQ( block(0,1), 8.);
    EXPECT_DOUBLE_EQ( block(0,2), 6.);
    EXPECT_DOUBLE_EQ( block(0,3), 15.);

    EXPECT_DOUBLE_EQ( block(2,0), 6.);
    EXPECT_DOUBLE_EQ( block(2,1), 9.);
    EXPECT_DOUBLE_EQ( block(2,2), 12.);
    EXPECT_DOUBLE_EQ( block(2,3), 18.);
  }
}

TEST(containers_multi_vector_eigen, mv_dot_mvEigen_ColMajorMatrix){
  Data d;
  // matrix with left layout
  using eig_mat_t = Eigen::Matrix<double, -1, -1, Eigen::ColMajor>;
  using mat_t = ::pressio::containers::Matrix<eig_mat_t>;
  doDot<mat_t>(d.A,d.B);
}

TEST(containers_multi_vector_eigen, mv_dot_mvEigen_RowMajorMatrix){
  Data d;
  // matrix with right layout
  using eig_mat_t = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
  using mat_t = ::pressio::containers::Matrix<eig_mat_t>;
  doDot<mat_t>(d.A,d.B);
}

TEST(containers_multi_vector_eigen, mv_dot_mvEigen_ColMajorMatrix_subspanStore){
  Data d;
  // matrix with left layout
  using eig_mat_t = Eigen::Matrix<double, -1, -1, Eigen::ColMajor>;
  using mat_t = ::pressio::containers::Matrix<eig_mat_t>;

  // make a matrix bit enough
  mat_t C(20, 25);
  C.setZero();

  // do operation on subspan of C from rows (10,13) and cols (15,19)
  doDot(d.A,d.B,C);

  // check orginal matrix is modified
  for (int i=0; i<20; i++){
    for (int j=0; j<25; j++){
      if (i<10) EXPECT_DOUBLE_EQ( C(i,j), 0.);
      if (i>12) EXPECT_DOUBLE_EQ( C(i,j), 0.);
      if (j<15) EXPECT_DOUBLE_EQ( C(i,j), 0.);
      if (j>18) EXPECT_DOUBLE_EQ( C(i,j), 0.);
    }
  }

  EXPECT_DOUBLE_EQ( C(10,15), 11.);
  EXPECT_DOUBLE_EQ( C(10,16), 8.);
  EXPECT_DOUBLE_EQ( C(10,17), 6.);
  EXPECT_DOUBLE_EQ( C(10,18), 15.);

  EXPECT_DOUBLE_EQ( C(12,15), 6.);
  EXPECT_DOUBLE_EQ( C(12,16), 9.);
  EXPECT_DOUBLE_EQ( C(12,17), 12.);
  EXPECT_DOUBLE_EQ( C(12,18), 18.);
}
