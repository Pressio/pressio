
#if not defined APPS_LSPG_UTILS_EPETRA_HPP_
#define APPS_LSPG_UTILS_EPETRA_HPP_

#include "CORE_ALL"
#include "SVD_BASIC"
#include "Epetra_MpiComm.h"

namespace rompp{ namespace apps{ namespace test{ namespace epetra{

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
      T nrows, T ncols,
      Epetra_MpiComm & Comm,
      const Epetra_Map & rowMap)
  -> rompp::core::MultiVector<Epetra_MultiVector>{

  rompp::core::MultiVector<Epetra_MultiVector> ADW(rowMap, ncols);
  // each process stores just its elements from A0
  int nMyElem = rowMap.NumMyElements();
  std::vector<int> myGel(nMyElem);
  rowMap.MyGlobalElements(myGel.data());
  for (int i=0; i<nMyElem; i++){
    int gi = myGel[i];
    for (int j=0; j<ncols; j++)
      ADW.replaceGlobalValue(gi, j, A0[gi][j]);
  }
  return ADW;
}

// template just to avoid having a cc file
template <typename T = int>
auto readBasis(
  std::string filename,
  T romSize, T numCell,
  Epetra_MpiComm & Comm,
  const Epetra_Map & rowMap)
  ->rompp::core::MultiVector<Epetra_MultiVector>
{
  std::vector<std::vector<double>> A0;
  readMatrixFromFile(filename, A0, romSize);
  // read basis into a MultiVector
  auto phi = convertFromVVecToMultiVec(A0, numCell, romSize, Comm, rowMap);
  //  phi.data()->Print(std::cout);
  return phi;
}

}}}}// end namespace rompp::apps::test::epetra

#endif
