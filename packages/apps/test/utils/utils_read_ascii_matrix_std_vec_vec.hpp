
#if not defined APPS_UTILS_READ_ASCII_MATRIX_STD_VEC_VEC_HPP_
#define APPS_UTILS_READ_ASCII_MATRIX_STD_VEC_VEC_HPP_

#include "CONTAINERS_ALL"
#include "SVD_BASIC"
#include <fstream>

namespace rompp{ namespace apps{ namespace test{ 

// template just to avoid having a cc file
template <typename T = int>
void readAsciiMatrixStdVecVec(std::string filename,
			std::vector<std::vector<double>> & A0,
			T ncols)
{
  assert( A0.empty() );
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

}}}// end namespace rompp::apps::test

#endif
