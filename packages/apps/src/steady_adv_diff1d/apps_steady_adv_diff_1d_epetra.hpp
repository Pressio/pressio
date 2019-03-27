
#ifndef ROM_APPS_STEADY_ADV_DIFF_1D_EPETRA_HPP_
#define ROM_APPS_STEADY_ADV_DIFF_1D_EPETRA_HPP_

#include "../../../CORE_ALL"
#include "Epetra_MpiComm.h"
#include <Epetra_config.h>
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "AztecOO_config.h"
#include "AztecOO.h"

namespace rompp{ namespace apps{

class SteadyAdvDiff1dEpetra{
protected:
  using nativeVec = Epetra_Vector;
  template<typename T> using rcp = std::shared_ptr<T>;
  using nativeMatrix	= Epetra_CrsMatrix;

  /* these types exposed because need to be detected */
public:
  using scalar_type	= double;
  using state_type	= Epetra_Vector;
  using residual_type	= state_type;

public:
  SteadyAdvDiff1dEpetra(Epetra_MpiComm & comm,
			std::vector<scalar_type> & mu,
			std::vector<scalar_type> & domain,
			std::vector<scalar_type> & bc1D)
    : comm_(comm), mu_(mu), domain_(domain), bc1D_(bc1D){}

  ~SteadyAdvDiff1dEpetra() = default;

  Epetra_Map const & getDataMap(){ return *contigMap_; };

public:
  void createMap();
  void setup();
  rcp<nativeMatrix> calculateLinearSystem();
  rcp<nativeVec> calculateForcingTerm();
  rcp<nativeVec> getStates();
  void calculateStates(rcp<nativeMatrix> A,
		       rcp<nativeVec> u,
		       rcp<nativeVec> f);
  rcp<nativeVec> calculateManufacturedForcing();
  void printStates(rcp<nativeVec> u);
  void compare2manufacturedStates(rcp<nativeVec> uapprox);
  double verifyImplementation(rcp<nativeMatrix> A);

protected:
  Epetra_MpiComm & comm_;
  std::vector<scalar_type> mu_;
  std::vector<scalar_type> domain_;
  std::vector<scalar_type> bc1D_;
  rcp<Epetra_Map> contigMap_;
  rcp<nativeMatrix> A_;
  int numGlobalNodes_;
  int *MyGlobalNodes_;
  int nodesPerProc_;
  double x_i_;
  rcp<nativeVec> u_;
  rcp<nativeVec> f_;
};

}} //namespace rompp::apps
#endif
