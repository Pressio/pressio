
#ifndef ROMPP_APPS_STEADY_LIN_ADV_DIFF_1D_EPETRA_HPP_
#define ROMPP_APPS_STEADY_LIN_ADV_DIFF_1D_EPETRA_HPP_

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

class SteadyLinAdvDiff1dEpetra{
protected:
  using nativeVec = Epetra_Vector;
  template<typename T> using rcp = std::shared_ptr<T>;
  using nativeMatrix  = Epetra_CrsMatrix;

public:
  /* these types exposed because need to be detected */
  using scalar_type = double;
  using state_type  = Epetra_Vector;
  using residual_type = state_type;
  using jacobian_type   = nativeMatrix;
public:
  SteadyLinAdvDiff1dEpetra(Epetra_MpiComm & comm,
      std::vector<scalar_type> & mu,
      std::vector<scalar_type> & domain,
      std::vector<scalar_type> & bc1D)
    : comm_(comm), mu_(mu), domain_(domain), bc1D_(bc1D),
      dx_{domain_[2]}, alpha_{mu_[0]}, beta_{mu_[1]},
      alphaOvDxSq_{alpha_/(dx_*dx_)},
      betaOvDx2_{beta_/(2.0*dx_)}
  {}

  ~SteadyLinAdvDiff1dEpetra() = default;

public:
  void createMap();
  Epetra_Map const & getDataMap(){ return *contigMap_; };
  void setup();
  void calculateLinearSystem() const;
  void calculateForcingTerm() const;
  int getNumGlobalNodes() const;
  rcp<nativeVec> getState() const;
  rcp<nativeVec> getGrid() const;
  rcp<nativeVec> getforcing() const;
  void solve();

  void printState() const{
    u_->Print( std::cout << std::setprecision(10) );
  }

public:
  void residual(const state_type & u,
    residual_type & rhs) const{
    /* compute jacobian and forcing term
     * (even though for this prob we do not need to
     * recompute every time, for sake of generality,
     * we keep it this way */
    calculateLinearSystem();
    calculateForcingTerm();

    A_->Multiply(false, u, rhs);
    // now, rhs = A*u so we just subtract f to obtain residual
    rhs.Update(-1., (*f_), 1.0);
  }

  residual_type residual(const state_type & u) const{
    Epetra_Vector R(*contigMap_);
    residual(u,R);
    return R;
  };

  // computes: C = Jac B where B is a multivector
  void applyJacobian(const state_type & y,
         const Epetra_MultiVector & B,
         Epetra_MultiVector & C) const
  {
    assert( A_->NumGlobalCols() == B.GlobalLength() );
    assert( C.GlobalLength() == A_->NumGlobalRows() );
    assert( C.NumVectors() == B.NumVectors() );
    /* compute jacobian (even though for this prob
     * we do not need to recompute every time,
     * for sake of generality, we keep it this way */
    calculateLinearSystem();
    A_->Multiply(false, B, C);
  }

//   // computes: A = Jac B where B is a multivector
  Epetra_MultiVector applyJacobian(const state_type & y,
             const Epetra_MultiVector & B) const{
    Epetra_MultiVector C( *contigMap_, B.NumVectors() );
    applyJacobian(y, B, C);
    return C;
  };

protected:
  Epetra_MpiComm & comm_;
  std::vector<scalar_type> mu_;
  std::vector<scalar_type> domain_;
  std::vector<scalar_type> bc1D_;
  scalar_type dx_{};
  scalar_type alpha_{};
  scalar_type beta_{};
  scalar_type alphaOvDxSq_{};
  scalar_type betaOvDx2_{};

  rcp<Epetra_Map> contigMap_;
  mutable rcp<nativeMatrix> A_;
  int numGlobalNodes_;
  int *MyGlobalNodes_;
  int nodesPerProc_;
  rcp<nativeVec> x_;
  mutable rcp<nativeVec> u_;
  mutable rcp<nativeVec> f_;
};

}} //namespace rompp::apps
#endif
