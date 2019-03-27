
#ifndef ROMPP_APPS_STEADY_ADV_DIFF_1D_EPETRA_HPP_
#define ROMPP_APPS_STEADY_ADV_DIFF_1D_EPETRA_HPP_

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
#include <cmath>

namespace rompp{ namespace apps{

class SteadyAdvDiff1dEpetra{
protected:
  using nativeVec = Epetra_Vector;
  template<typename T> using rcp = std::shared_ptr<T>;
  using nativeMatrix	= Epetra_CrsMatrix;

public:
  /* these types exposed because need to be detected */
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

public:
  void createMap();
  Epetra_Map const & getDataMap(){ return *contigMap_; };
  void setup();
  void calculateLinearSystem();
  void calculateForcingTerm();
  int getNumGlobalNodes() const;
  rcp<nativeVec> getState() const;
  rcp<nativeVec> getGrid() const;
  void solve();
  void printState() const;

public:
  void residual(const state_type & u,
		residual_type & rhs) const{
    // compute residual and store into rhs
  }

  residual_type residual(const state_type & u) const{
    /* this should create a vector, compure residual and return it */

    Epetra_Vector R(*contigMap_);
    // residual(u,R,t);
    return R;
  }

  // computes: A = Jac B where B is a multivector
  void applyJacobian(const state_type & y,
		     const Epetra_MultiVector & B,
		     Epetra_MultiVector & A) const{
    // assert( Jac_->NumGlobalCols() == B.GlobalLength() );
    // assert( A.GlobalLength() == Jac_->NumGlobalRows() );
    // assert( A.NumVectors() == B.NumVectors() );
    // // compute jacobian
    // jacobian(y, *Jac_, t);
    // Jac_->Multiply(false, B, A);
  }

  // computes: A = Jac B where B is a multivector
  Epetra_MultiVector applyJacobian(const state_type & y,
  				   const Epetra_MultiVector & B) const{
    Epetra_MultiVector C( *contigMap_, B.NumVectors() );
    // applyJacobian(y, B, C, t);
    return C;
  }

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
  rcp<nativeVec> x_;
  rcp<nativeVec> u_;
  rcp<nativeVec> f_;
};

}} //namespace rompp::apps
#endif
