
#if not defined APPS_UTILS_EIGEN_HPP_
#define APPS_UTILS_EIGEN_HPP_

#include "CONTAINERS_ALL"
#include "SVD_BASIC"
#include "utils_read_ascii_matrix_std_vec_vec.hpp"

namespace rompp{ namespace apps{ namespace test{ namespace eigen{

// template just to avoid having a cc file
template <typename T = int>
auto convertFromVVecToMultiVec(
      const std::vector<std::vector<double>> & A0,
      T nrows, T ncols)
  -> rompp::containers::MultiVector<Eigen::MatrixXd>{

  rompp::containers::MultiVector<Eigen::MatrixXd> ADW(nrows, ncols);

  for (int i=0; i<nrows; i++){
    for (int j=0; j<ncols; j++)
      ADW(i,j) = A0[i][j];
  }
  return ADW;
}

// template just to avoid having a cc file
template <typename T = int>
auto readBasis(
  std::string filename,
  T romSize, T numCell)
  ->rompp::containers::MultiVector<Eigen::MatrixXd>
{
  std::vector<std::vector<double>> A0;
  ::rompp::apps::test::readAsciiMatrixStdVecVec(filename, A0, romSize);
  // read basis into a MultiVector
  auto phi = convertFromVVecToMultiVec(A0, numCell, romSize);
  //  phi.data()->Print(std::cout);
  return phi;
}

}}}}// end namespace rompp::apps::test::eigen

#endif
