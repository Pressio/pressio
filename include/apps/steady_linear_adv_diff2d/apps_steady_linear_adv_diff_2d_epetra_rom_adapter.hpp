/*
//@HEADER
// ************************************************************************
//
// apps_steady_linear_adv_diff_2d_epetra_rom_adapter.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
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

#ifndef APPS_STEADY_LINEAR_ADV_DIFF2D_APPS_STEADY_LINEAR_ADV_DIFF_2D_EPETRA_ROM_ADAPTER_HPP_
#define APPS_STEADY_LINEAR_ADV_DIFF2D_APPS_STEADY_LINEAR_ADV_DIFF_2D_EPETRA_ROM_ADAPTER_HPP_

#include "apps_steady_linear_adv_diff_2d_epetra.hpp"

namespace pressio{ namespace apps{

class SteadyLinAdvDiff2dEpetraRomAdapter
{
public:
  /* these types exposed because need to be detected */
  using scalar_type	= double;
  using state_type	= Epetra_Vector;
  using residual_type  = Epetra_Vector;
  using dense_matrix_type = Epetra_MultiVector;

private:
  using velocity_type = Epetra_Vector;

public:
  template <typename ... Args>
  SteadyLinAdvDiff2dEpetraRomAdapter(Args&& ... args)
    : appObj_{std::forward<Args>(args)...}
  {}

public:
  Epetra_Map const & getDataMap()const {
    return appObj_.getDataMap();
  };

  void printStateToFile(std::string fileName,
			state_type & T){
    auto x = appObj_.getX();
    auto y = appObj_.getY();
    auto dofPerProc = appObj_.getNumLocalDofs();

    std::ofstream file;
    file.open( fileName );
    for(auto i=0; i < dofPerProc; i++){
      file << std::fixed << std::setprecision(14) <<
	(*x)[i] << " " << (*y)[i] << " " << T[i];
      file << std::endl;
    }
    file.close();
  }

  std::shared_ptr<state_type>
  getState() const {
    return appObj_.getState();
  }

  residual_type createResidual() const{
    residual_type R( appObj_.getDataMap() );
    return R;
  };

  // computes: A = Jac B where B is a multivector
  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const
  {
    dense_matrix_type C( appObj_.getDataMap(), B.NumVectors() );
    return C;
  };

  void residual(const state_type & yState, residual_type & rhs) const
  {
    appObj_.assembleMatrix();
    appObj_.fillRhs();
    auto A = appObj_.getMatrix();
    A->Multiply(false, yState, rhs);
    // now, rhs = A*u, subtract forcing term because forcing
    // was calculated by moving the terms to the rhs
    auto f = appObj_.getForcing();
    rhs.Update(-1., (*f), 1.0);
  }

  // computes: C = Jac B where B is a multivector
  void applyJacobian(const state_type & yState,
		     const dense_matrix_type & B,
		     dense_matrix_type & C) const
  {
    appObj_.assembleMatrix();
    auto A = appObj_.getMatrix();
    assert( A->NumGlobalCols() == B.GlobalLength() );
    assert( C.GlobalLength() == A->NumGlobalRows() );
    assert( C.NumVectors() == B.NumVectors() );
    A->Multiply(false, B, C);
  }

private:
  SteadyLinAdvDiff2dEpetra appObj_;

};

}} //namespace pressio::apps
#endif  // APPS_STEADY_LINEAR_ADV_DIFF2D_APPS_STEADY_LINEAR_ADV_DIFF_2D_EPETRA_ROM_ADAPTER_HPP_
