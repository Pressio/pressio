
#ifndef ROM_TEST_LAPLACE1D_EPETRA_HPP_
#define ROM_TEST_LAPLACE1D_EPETRA_HPP_

#include "CORE_ALL"
#include "ODE_ALL"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "AztecOO_config.h"
#include "AztecOO.h"

class Laplace1dEpetra{
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
  Laplace1dEpetra(Epetra_MpiComm * comm, 
		  std::vector<scalar_type> mu,
		  std::vector<scalar_type> domain, 
		  std::vector<scalar_type> bc1D)
    : comm_(comm), mu_(mu), domain_(domain), bc1D_(bc1D){}

  ~Laplace1dEpetra() = default;
  
  Epetra_Map const & getDataMap(){ return *contigMap; };

public:
  rcp<nativeMatrix> calculateLinearSystem();
  rcp<nativeVec> calculateForcingTerm();
  rcp<nativeVec> createStates();
  void solveForStates(rcp<nativeMatrix> A, rcp<nativeVec> u, rcp<nativeVec> f);
  rcp<nativeVec> calcManufacturedForcing();
  void printStates(rcp<nativeVec> u);
  void compare2manufacturedStates(rcp<nativeVec> uapprox);
  double verifyImplementation(rcp<nativeMatrix> A);

protected:
  Epetra_MpiComm * comm_;
  rcp<Epetra_Map> contigMap;
  std::vector<scalar_type> domain_;
  std::vector<scalar_type> mu_;
  std::vector<scalar_type> bc1D_;
  rcp<nativeMatrix> A;
  int numGlobalNodes;
  int *MyGlobalNodes;
  int nodesPerProc; 
  int *NumNz;
  double *Values;
  int *Indices;
  double x_i;
  double two;
  double one;
  int NumEntries;
  rcp<nativeVec> state_u;
  rcp<nativeVec> f;
};

#endif
