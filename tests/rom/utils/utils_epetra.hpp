
#ifndef ROM_TEST_UTILS_EPETRA_HPP_
#define ROM_TEST_UTILS_EPETRA_HPP_

#include "pressio_utils.hpp"
#include "Epetra_MpiComm.h"

namespace pressio{ namespace rom{ namespace test{ namespace epetra{

// template just to avoid having a cc file
template <typename T = int>
auto convertFromVVecToMultiVec(
      const std::vector<std::vector<double>> & A0,
      const T nrows, const T ncols,
      const Epetra_MpiComm & Comm,
      const Epetra_Map & rowMap)
  -> pressio::containers::MultiVector<Epetra_MultiVector>{

  pressio::containers::MultiVector<Epetra_MultiVector> ADW(rowMap, ncols);
  // each process stores just its elements from A0
  int nMyElem = rowMap.NumMyElements();
  std::vector<int> myGel(nMyElem);
  rowMap.MyGlobalElements(myGel.data());
  for (int i=0; i<nMyElem; i++){
    int gi = myGel[i];
    for (int j=0; j<ncols; j++)
      ADW.data()->ReplaceGlobalValue(gi, j, A0[gi][j]);
  }
  return ADW;
}

// template just to avoid having a cc file
template <typename T = int>
auto readBasis(
  std::string filename,
  T romSize, T numCell,
  const Epetra_MpiComm & Comm,
  const Epetra_Map & rowMap)
  ->pressio::containers::MultiVector<Epetra_MultiVector>
{
  std::vector<std::vector<double>> A0;
  ::pressio::utils::readAsciiMatrixStdVecVec(filename, A0, romSize);
  // read basis into a MultiVector
  auto phi = convertFromVVecToMultiVec(A0, numCell, romSize, Comm, rowMap);
  //  phi.data()->Print(std::cout);
  return phi;
}

}}}}// end namespace pressio::rom::test::epetra

#endif

