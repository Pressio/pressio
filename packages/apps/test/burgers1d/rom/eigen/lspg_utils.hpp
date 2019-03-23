
#if not defined APPS_LSPG_UTILS_EIGEN_HPP_
#define APPS_LSPG_UTILS_EIGEN_HPP_

#include "CORE_ALL"
#include "SVD_BASIC"

namespace rompp{ namespace apps{ namespace test{ namespace eigen{

// template just to avoid having a cc file
template <typename T = int>
void readMatrixFromFile(std::string filename,
			std::vector<std::vector<double>> & A0,
			T ncols){

  std::ifstream source;
  source.open( filename, std::ios_base::in);
  std::string line, colv;
  std::vector<double> tmpv(ncols);
  while (std::getline(source, line) ){
    std::istringstream in(line);
    for (int i=0; i<ncols; i++){
      in >> colv;
      tmpv[i] = atof(colv.c_str());
    }
    A0.emplace_back(tmpv);
  }
  source.close();
}

// template just to avoid having a cc file
template <typename T = int>
auto convertFromVVecToMultiVec(
      const std::vector<std::vector<double>> & A0,
      T nrows, T ncols)
  -> rompp::core::MultiVector<Eigen::MatrixXd>{

  rompp::core::MultiVector<Eigen::MatrixXd> ADW(nrows, ncols);

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
  ->rompp::core::MultiVector<Eigen::MatrixXd>
{
  std::vector<std::vector<double>> A0;
  readMatrixFromFile(filename, A0, romSize);
  // read basis into a MultiVector
  auto phi = convertFromVVecToMultiVec(A0, numCell, romSize);
  //  phi.data()->Print(std::cout);
  return phi;
}

}}}}// end namespace rompp::apps::test::eigen

#endif
