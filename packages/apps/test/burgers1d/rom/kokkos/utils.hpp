
#if not defined APPS_LSPG_UTILS_KOKKOS_HPP_
#define APPS_LSPG_UTILS_KOKKOS_HPP_

#include "UTILS_ALL"
#include "CONTAINERS_ALL"

namespace pressio{ namespace apps{ namespace test{ namespace kokkos{

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

template <typename T = int>
auto convertFromVVecToMultiVec(
      const std::vector<std::vector<double>> & A0,
      T nrows, T ncols)
  -> pressio::containers::MultiVector<Eigen::MatrixXd>{

  pressio::containers::MultiVector<Eigen::MatrixXd> ADW(nrows, ncols);

  for (int i=0; i<nrows; i++){
    for (int j=0; j<ncols; j++)
      ADW(i,j) = A0[i][j];
  }
  return ADW;
}

template <typename phi_d_t>
void readBasis(std::string filename,
	       int romSize,
	       int numCell,
	       phi_d_t phi_d)
{
  std::vector<std::vector<double>> A0;
  readMatrixFromFile(filename, A0, romSize);

  using view_h_t = typename phi_d_t::HostMirror;
  view_h_t phi_h("phi_h", numCell, romSize);

  for (int i=0; i<numCell; i++){
    for (int j=0; j<romSize; j++)
      phi_h(i,j) = A0[i][j];
  }
  // deep copy to device
  Kokkos::deep_copy(phi_d, phi_h);
}

}}}}// end namespace pressio::apps::test::eigen

#endif
