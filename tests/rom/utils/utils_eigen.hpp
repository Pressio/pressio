
#if not defined ROM_TEST_UTILS_EIGEN_HPP_
#define ROM_TEST_UTILS_EIGEN_HPP_

#include "pressio_utils.hpp"

namespace pressio{ namespace rom{ namespace test{ namespace eigen{

template <typename T = int>
auto convertFromVVecToMultiVec(
      const std::vector<std::vector<double>> & A0,
      T nrows, T ncols)
  -> pressio::containers::MultiVector<Eigen::MatrixXd>{

  using eig_mat = Eigen::MatrixXd;
  pressio::containers::MultiVector<eig_mat> ADW(nrows, ncols);

  for (int i=0; i<nrows; i++){
    for (int j=0; j<ncols; j++)
      ADW(i,j) = A0[i][j];
  }
  return ADW;
}

template <typename T = int>
auto readBasis(
  std::string filename,
  T romSize, T numCell)
  ->pressio::containers::MultiVector<Eigen::MatrixXd>
{
  std::vector<std::vector<double>> A0;
  ::pressio::utils::readAsciiMatrixStdVecVec(filename, A0, romSize);
  // read basis into a MultiVector
  auto phi = convertFromVVecToMultiVec(A0, numCell, romSize);
  //  phi.data()->Print(std::cout);
  return phi;
}

}}}}// end namespace pressio::rom::test::eigen

#endif
