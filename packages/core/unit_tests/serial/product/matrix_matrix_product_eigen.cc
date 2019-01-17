
#include <gtest/gtest.h>
#include "CORE_ALL"

TEST(core_matrix_matrix_product, eigenDenseDense){
  using namespace rompp;

  using nat_t = Eigen::MatrixXd;
  using myA_t = core::Matrix<nat_t>;
  nat_t a(3,4);
  a << 1.,2.,3.,4., 4.,3.,2.,1., 1.,2.,3.,4.;
  myA_t A(a);
  //  std::cout << *A.data() << "\n";

  using nat2_t = Eigen::MatrixXd;
  using myB_t = core::Matrix<nat2_t>;
  nat2_t b(4,2); b << 1.,2.,3.,3.,4.,4.,2.,2.;
  myB_t B(b);
  //  std::cout << *B.data();

  static_assert(
   core::meta::is_dense_matrix_wrapper_eigen<myA_t>::value,"");
  
  auto C = core::ops::product(A,B);
  EXPECT_EQ( C.cols(), B.cols());
  EXPECT_EQ( C.rows(), A.rows());
  EXPECT_DOUBLE_EQ( C(0,0), 27.0);
  EXPECT_DOUBLE_EQ( C(1,0), 23.0);

  core::Matrix<Eigen::Matrix<double,3,2>> C2;
  core::ops::product(A,B,C2);
  EXPECT_DOUBLE_EQ( C2(0,0), 27.0);
  EXPECT_DOUBLE_EQ( C2(1,0), 23.0);

}//end TEST


TEST(core_matrix_matrix_product, eigenSparseSparse){  
  using namespace rompp;  
  // create sparse matrix
  using mymat_t = core::Matrix<
    Eigen::SparseMatrix<double,Eigen::RowMajor,int>>;
  mymat_t A(4,3);
  {
  using veci = std::vector<int>;
  using vecd = std::vector<double>;
  std::map<int,std::pair<vecd,veci>> DD;
  DD[0] = std::make_pair<vecd,veci>({1.,2.}, {0,2});
  DD[1] = std::make_pair<vecd,veci>({2.,3.}, {0,2});
  DD[2] = std::make_pair<vecd,veci>({1}, {2});
  DD[3] = std::make_pair<vecd,veci>({1}, {1});
  for (auto const & it : DD){
    int ind = it.first;
    vecd data = std::get<0>(it.second);
    veci cols = std::get<1>(it.second);
    A.insertValues(ind,(int)data.size(),data.data(),cols.data());
  }}
  
  // create sparse matrix
  using mymat_t = core::Matrix<
    Eigen::SparseMatrix<double,Eigen::RowMajor,int>>;
  mymat_t B(3,2);
  {
  using veci = std::vector<int>;
  using vecd = std::vector<double>;
  std::map<int,std::pair<vecd,veci>> DD;
  DD[0] = std::make_pair<vecd,veci>({1}, {0});
  DD[1] = std::make_pair<vecd,veci>({1.},{1});
  for (auto const & it : DD){
    int ind = it.first;
    vecd data = std::get<0>(it.second);
    veci cols = std::get<1>(it.second);
    B.insertValues(ind,(int)data.size(),data.data(),cols.data());
  }}

  auto C = core::ops::product(A,B);
  EXPECT_DOUBLE_EQ( C(0,0), 1.);
  EXPECT_DOUBLE_EQ( C(0,1), 0.0);
  EXPECT_DOUBLE_EQ( C(1,0), 2.0);
  EXPECT_DOUBLE_EQ( C(1,1), 0.0);
  EXPECT_DOUBLE_EQ( C(2,0), 0.0);
  EXPECT_DOUBLE_EQ( C(2,1), 0.0);
  EXPECT_DOUBLE_EQ( C(3,0), 0.0);
  EXPECT_DOUBLE_EQ( C(3,1), 1.0);

  mymat_t C2(4,2);
  core::ops::product(A,B,C2);
  EXPECT_DOUBLE_EQ( C2(0,0), 1.);
  EXPECT_DOUBLE_EQ( C2(0,1), 0.0);
  EXPECT_DOUBLE_EQ( C2(1,0), 2.0);
  EXPECT_DOUBLE_EQ( C2(1,1), 0.0);
  EXPECT_DOUBLE_EQ( C2(2,0), 0.0);
  EXPECT_DOUBLE_EQ( C2(2,1), 0.0);
  EXPECT_DOUBLE_EQ( C2(3,0), 0.0);
  EXPECT_DOUBLE_EQ( C2(3,1), 1.0);
  
}//end TEST




TEST(core_matrix_matrix_product, eigenSparseDense){
  using namespace rompp;  
  // create sparse matrix
  using mymat_t = core::Matrix<
    Eigen::SparseMatrix<double,Eigen::RowMajor,int>>;
  mymat_t A(4,3);
  {
  using veci = std::vector<int>;
  using vecd = std::vector<double>;
  std::map<int,std::pair<vecd,veci>> DD;
  DD[0] = std::make_pair<vecd,veci>({1.,2.}, {0,2});
  DD[1] = std::make_pair<vecd,veci>({2.,3.}, {0,2});
  DD[2] = std::make_pair<vecd,veci>({1}, {2});
  DD[3] = std::make_pair<vecd,veci>({1}, {1});
  for (auto const & it : DD){
    int ind = it.first;
    vecd data = std::get<0>(it.second);
    veci cols = std::get<1>(it.second);
    A.insertValues(ind,(int)data.size(),data.data(),cols.data());
  }}
  
  // create dense matrix
  using nat_t = Eigen::MatrixXd;
  using myB_t = core::Matrix<nat_t>;
  nat_t b(3,4);
  b << 1.,2.,3.,4., 4.,3.,2.,1., 1.,2.,3.,4.;
  myB_t B(b);
  //  std::cout << *A.data() << "\n";

  // do product  
  auto C = core::ops::product(A,B);
  //  std::cout << *C.data();
  EXPECT_DOUBLE_EQ( C(0,0), 3.);
  EXPECT_DOUBLE_EQ( C(0,1), 6.);
  EXPECT_DOUBLE_EQ( C(0,2), 9.);
  EXPECT_DOUBLE_EQ( C(0,3), 12.);

  EXPECT_DOUBLE_EQ( C(1,0), 5.);
  EXPECT_DOUBLE_EQ( C(1,1), 10.);
  EXPECT_DOUBLE_EQ( C(1,2), 15.);
  EXPECT_DOUBLE_EQ( C(1,3), 20.);

  EXPECT_DOUBLE_EQ( C(2,0), 1.);
  EXPECT_DOUBLE_EQ( C(2,1), 2.);
  EXPECT_DOUBLE_EQ( C(2,2), 3.);
  EXPECT_DOUBLE_EQ( C(2,3), 4.);

  EXPECT_DOUBLE_EQ( C(3,0), 4.);
  EXPECT_DOUBLE_EQ( C(3,1), 3.);
  EXPECT_DOUBLE_EQ( C(3,2), 2.);
  EXPECT_DOUBLE_EQ( C(3,3), 1.);

  core::Matrix<Eigen::MatrixXd> C2(4,4);
  core::ops::product(A,B,C2);
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      EXPECT_DOUBLE_EQ( C2(i,j), C(i,j) );
  
}//end TEST
