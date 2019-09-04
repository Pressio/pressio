/*
//@HEADER
// ************************************************************************
//
// apps_unsteady_linear_adv_diff_1d_epetra.hpp
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

#ifndef PRESSIO_APPS_UNSTEADY_LINEAR_ADV_DIFF_1D_EPETRA_HPP_
#define PRESSIO_APPS_UNSTEADY_LINEAR_ADV_DIFF_1D_EPETRA_HPP_

#include "../apps_ConfigDefs.hpp"

#ifdef HAVE_TRILINOS
#include "../../../CONTAINERS_ALL"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include <cmath>
#include "../steady_linear_adv_diff1d/apps_steady_linear_adv_diff_1d_epetra.hpp"

namespace pressio{ namespace apps{
class UnsteadyLinAdvDiff1dEpetra: public SteadyLinAdvDiff1dEpetra{
protected:
  using nativeVec = Epetra_Vector;
  using nativeMatrix  = Epetra_CrsMatrix;

public:
  UnsteadyLinAdvDiff1dEpetra(const Epetra_MpiComm & comm,
			     const std::vector<scalar_type> & mu,
			     const std::vector<scalar_type> & domain,
			     const std::vector<scalar_type> & bc1D)
    : SteadyLinAdvDiff1dEpetra(comm, mu, domain, bc1D){}
  ~UnsteadyLinAdvDiff1dEpetra() = default;

public:
  void unsteadySetup();

  rcp<nativeVec> getInitialState() const;

  void velocity(const state_type & u,
		const scalar_type /* t*/,
    velocity_type & rhs) const;

  velocity_type velocity(const state_type & u,
			 const scalar_type t) const{
    Epetra_Vector R(*contigMap_);
    velocity(u, t, R);
    return R;
  }

  void applyJacobian(const state_type & y,
		     const Epetra_MultiVector & B,
		     scalar_type /*t*/,
         Epetra_MultiVector &A) const{
    SteadyLinAdvDiff1dEpetra::applyJacobian(y, B, A);
    A.Scale(-1.0);
  }

  Epetra_MultiVector applyJacobian(const state_type &y,
				   const Epetra_MultiVector &B,
				   scalar_type t) const{
    Epetra_MultiVector C( *contigMap_, B.NumVectors());
    applyJacobian(y, B, t, C);
    return C;
  }

protected:
  mutable rcp<nativeVec> U0_;
};

}}
#endif
#endif
