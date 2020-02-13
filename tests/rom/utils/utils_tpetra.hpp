
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#if not defined ROM_TEST_UTILS_TPETRA_HPP_
#define ROM_TEST_UTILS_TPETRA_HPP_

#include "pressio_utils.hpp"
#include <Tpetra_Core.hpp>

namespace pressio{ namespace rom{ namespace test{ namespace tpetra{

template <typename rcpcomm_t, typename rcpmap_t>
auto convertFromVVecToMultiVec(const std::vector<std::vector<double>> & A0,
			       int nrows, int ncols,
			       rcpcomm_t Comm,
			       rcpmap_t rowMap)
  -> pressio::containers::MultiVector<Tpetra::MultiVector<>>{

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

  pressio::containers::MultiVector<Tpetra::MultiVector<>> ADW(AD);
  return ADW;
}


template <typename comm_t, typename map_t>
auto readBasis(
  std::string filename,
  int romSize, int numCell,
  comm_t Comm,
  const map_t rowMap)
  ->pressio::containers::MultiVector<Tpetra::MultiVector<>>
{
  std::vector<std::vector<double>> A0;
  ::pressio::utils::readAsciiMatrixStdVecVec(filename, A0, romSize);
  // read basis into a MultiVector
  return convertFromVVecToMultiVec(A0, numCell, romSize, Comm, rowMap);
}

}}}}// end namespace pressio::rom::test::tpetra

#endif
#endif

