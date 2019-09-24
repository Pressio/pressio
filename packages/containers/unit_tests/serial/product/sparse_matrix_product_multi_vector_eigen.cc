
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"

TEST(containers_ops_eigen, sparseMatProdMultivector){
  using namespace pressio;
  
  using esp_t = Eigen::SparseMatrix<double,Eigen::RowMajor,int>;
  using mymat_t = containers::Matrix<esp_t>;
    
  // create sparse matrix
  mymat_t A(6,4);
  {
  using veci = std::vector<int>;
  using vecd = std::vector<double>;
  std::map<int,std::pair<vecd,veci>> DD;
  DD[0] = std::make_pair<vecd,veci>({1.,2.,3.}, {0,1,2});
  DD[1] = std::make_pair<vecd,veci>({3.,2.,1.}, {0,1,2});
  DD[2] = std::make_pair<vecd,veci>({1}, {2});
  DD[3] = std::make_pair<vecd,veci>({1.,1.}, {1,3});
  DD[4] = std::make_pair<vecd,veci>({1}, {0});
  DD[5] = std::make_pair<vecd,veci>({1.,1.}, {1,2});
  for (auto const & it : DD){
    int ind = it.first;
    vecd data = std::get<0>(it.second);
    veci cols = std::get<1>(it.second);
    A.insertValues(ind,(int)data.size(),
		   data.data(),
		   cols.data());
  }}

  //construct multivector
  using myMV_t = containers::MultiVector<Eigen::MatrixXd>;
  myMV_t mv(4,3);
  mv(0,0) = 1.; mv(0,1) = 2.; mv(0,2) = 3.;
  mv(1,0) = 3.; mv(1,1) = 2.; mv(1,2) = 1.;
  mv(2,0) = 0.; mv(2,1) = 0.; mv(2,2) = 1.;
  mv(3,0) = 0.; mv(3,1) = 1.; mv(3,2) = 0.;
  
  auto c1 = pressio::containers::ops::product<mymat_t, myMV_t, myMV_t>(A,mv);
  EXPECT_DOUBLE_EQ( c1(0,0), 7.);
  EXPECT_DOUBLE_EQ( c1(0,1), 6.);
  EXPECT_DOUBLE_EQ( c1(0,2), 8.);

  EXPECT_DOUBLE_EQ( c1(1,0), 9.);
  EXPECT_DOUBLE_EQ( c1(1,1), 10.);
  EXPECT_DOUBLE_EQ( c1(1,2), 12.);

  EXPECT_DOUBLE_EQ( c1(2,0), 0.);
  EXPECT_DOUBLE_EQ( c1(2,1), 0.);
  EXPECT_DOUBLE_EQ( c1(2,2), 1.);

  EXPECT_DOUBLE_EQ( c1(3,0), 3.);
  EXPECT_DOUBLE_EQ( c1(3,1), 3.);
  EXPECT_DOUBLE_EQ( c1(3,2), 1.);

  EXPECT_DOUBLE_EQ( c1(4,0), 1.);
  EXPECT_DOUBLE_EQ( c1(4,1), 2.);
  EXPECT_DOUBLE_EQ( c1(4,2), 3.);

  EXPECT_DOUBLE_EQ( c1(5,0), 3.);
  EXPECT_DOUBLE_EQ( c1(5,1), 2.);
  EXPECT_DOUBLE_EQ( c1(5,2), 2.);
  
  pressio::containers::Matrix<Eigen::MatrixXd> c(6,3);
  pressio::containers::ops::product(A,mv,c);
  for (auto i=0; i<c1.data()->rows(); i++)
    for (auto j=0; j<c1.data()->cols(); j++)
      EXPECT_DOUBLE_EQ( c(i,j), c1(i,j));
}
