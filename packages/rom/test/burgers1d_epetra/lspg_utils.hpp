
#if not defined LSPG_UTILS_HPP_
#define LSPG_UTILS_HPP_

#include "CORE_ALL"
#include "SVD_BASIC"
#include "Epetra_MpiComm.h"


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


// template just to avoid having a cc file
template <typename T = int>
auto convertFromVVecToMultiVec(const std::vector<std::vector<double>> & A0,
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
auto getBasis(T romSize, T numCell,
	 Epetra_MpiComm & Comm, const Epetra_Map & rowMap)
  ->rompp::core::MultiVector<Epetra_MultiVector>
{
  // //-----------------------------------------
  // // run forward Euler and collect snapshots
  //  int numSnaps = 20;
  // using coll_t = snaps_collector<app_state_w_t>;
  // coll_t obsO(numSnaps, appobj.getDataMap(), y0n);
  // {
  //   // types used for ode = wrappers types
  //   app_state_w_t y(y0n);
  //   app_res_w_t r(r0n);
  //   //stepper
  //   using stepper_t = rompp::ode::ExplicitStepper<
  //     rompp::ode::ExplicitSteppersEnum::Euler, app_state_w_t, app_t>;
  //   stepper_t stepperObj(appobj, y, r);
  //   // integrate in time
  //   rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, nSteps, obsO);
  //   //y.data()->Print(std::cout << std::setprecision(13));
  //   //    obsO.printAll();
  // }

  /// create svd solver
  int numBasis = romSize;
  // using snaps_storage_t = typename coll_t::snaps_storage_t;
  // rompp::svd::Solver<snaps_storage_t,
  // 		     rompp::core::MultiVector,
  // 		     rompp::core::MultiVector,
  // 		     std::vector<double>> svdO;
  // svdO.compute<rompp::svd::svdType::truncated>(obsO.snapshots_, numBasis);
  // using phi_t = rompp::core::MultiVector<Epetra_MultiVector>;
  // const phi_t & phi = svdO.cRefLeftSingularVectors();
  // phi.data()->Print(std::cout << std::setprecision(13));

  std::vector<std::vector<double>> A0;
  //  using phi_t = rompp::core::MultiVector<Epetra_MultiVector>;
  auto fin = "U" + std::to_string(numBasis) +
    "_ncell" + std::to_string(numCell) + "_t35_dt001.txt";
  readMatrixFromFile(fin, A0, numBasis);
  // read basis into a MultiVector
  auto phi = convertFromVVecToMultiVec(A0, numCell, numBasis, Comm, rowMap);
  //  phi.data()->Print(std::cout);
  return phi;
}

}}}// end namespace rompp::rom::test

#endif
