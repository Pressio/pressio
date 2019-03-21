
#if not defined LSPG_UTILS_HPP_
#define LSPG_UTILS_HPP_

#include "CORE_ALL"
#include "SVD_BASIC"

namespace rompp{ namespace rom{ namespace test{

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


template <typename rcpcomm_t, typename rcpmap_t>
auto convertFromVVecToMultiVec(const std::vector<std::vector<double>> & A0,
			       int nrows, int ncols,
			       rcpcomm_t Comm,
			       rcpmap_t rowMap)
  -> rompp::core::MultiVector<Tpetra::MultiVector<>>{

  Tpetra::MultiVector<> AD(rowMap, ncols);

  auto nMyElem = rowMap->getNodeNumElements();

  std::vector<int> myGel(nMyElem);
  auto minGId = rowMap->getMinGlobalIndex();
  myGel.resize(nMyElem);
  std::iota(myGel.begin(), myGel.end(), minGId);

  for (decltype(nMyElem) i=0; i<nMyElem; i++){
    int gi = myGel[i];
    for (int j=0; j<ncols; j++)
      AD.replaceGlobalValue(gi, j, A0[gi][j]);
  }

  rompp::core::MultiVector<Tpetra::MultiVector<>> ADW(AD);
  return ADW;
}


template <typename comm_t, typename map_t>
auto getBasis(int romSize, int numCell,
	      comm_t Comm,
	      const map_t rowMap)
  ->rompp::core::MultiVector<Tpetra::MultiVector<>>
{
  int numBasis = romSize;
  std::vector<std::vector<double>> A0;

  auto fin = "U" + std::to_string(numBasis) +
    "_ncell" + std::to_string(numCell) + "_t35_dt001.txt";

  readMatrixFromFile(fin, A0, numBasis);
  return convertFromVVecToMultiVec(A0, numCell, numBasis, Comm, rowMap);
}

}}}// end namespace rompp::rom::test

#endif
