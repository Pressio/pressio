/*
//@HEADER
// ************************************************************************
//
// apps_steady_linear_adv_diff_1d_epetra.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the 
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions 
// are met:
//
// 1. Redistributions of source code must retain the above copyright 
// notice, this list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its 
// contributors may be used to endorse or promote products derived 
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef PRESSIO_APPS_STEADY_LIN_ADV_DIFF_1D_EPETRA_HPP_
#define PRESSIO_APPS_STEADY_LIN_ADV_DIFF_1D_EPETRA_HPP_

#include "../apps_ConfigDefs.hpp"

#ifdef HAVE_TRILINOS

#include "../../../CONTAINERS_ALL"
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

namespace pressio{ namespace apps{

class SteadyLinAdvDiff1dEpetra{
protected:
  using nativeVec = Epetra_Vector;
  template<typename T> using rcp = std::shared_ptr<T>;
  using nativeMatrix  = Epetra_CrsMatrix;

public:
  /* these types exposed because need to be detected */
  using scalar_type = double;
  using state_type  = Epetra_Vector;
  using velocity_type = state_type;
  using jacobian_type   = nativeMatrix;

public:
  SteadyLinAdvDiff1dEpetra(const Epetra_MpiComm & comm,
			   const std::vector<scalar_type> & mu,
			   const std::vector<scalar_type> & domain,
			   const std::vector<scalar_type> & bc1D)
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
  rcp<nativeVec> getRHSforce() const;
  rcp<nativeMatrix> getLHSmatrix() const;
  void solve();

  void printState() const{
    u_->Print( std::cout << std::setprecision(10) );
  }

public:
  void velocity(const state_type & u,
    velocity_type & rhs) const{
    /* compute jacobian and forcing term
     * (even though for this prob we do not need to
     * recompute every time, for sake of generality,
     * we keep it this way */
    calculateLinearSystem();
    calculateForcingTerm();

    A_->Multiply(false, u, rhs);
    // now, rhs = A*u so we just subtract f to obtain velocity
    rhs.Update(-1., (*f_), 1.0);
  }

  velocity_type velocity(const state_type & u) const{
    Epetra_Vector R(*contigMap_);
    velocity(u,R);
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
  const Epetra_MpiComm & comm_;
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

}} //namespace pressio::apps
#endif
#endif
