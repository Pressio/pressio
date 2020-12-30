/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_epetra_preconditioned.hpp
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

#ifndef APPS_BURGERS1D_APPS_BURGERS1D_EPETRA_PRECONDITIONED_HPP_
#define APPS_BURGERS1D_APPS_BURGERS1D_EPETRA_PRECONDITIONED_HPP_

#include "apps_burgers1d_epetra.hpp"

namespace pressio{ namespace apps{

class EpetraIdentityPreconditioner
{
  using scalar_type = typename Burgers1dEpetra::scalar_type;
  using state_type  = typename Burgers1dEpetra::state_type;
  using velocity_type = typename Burgers1dEpetra::velocity_type;

public:
  void applyPreconditioner(const state_type & yState,
            const scalar_type & time,
            velocity_type & rhs) const 
  {
    // do nothing, preconditioner is identity
    std::cout << "identiy precond" << std::endl;
  }

  void applyPreconditioner(const state_type & yState,
            const scalar_type & time,
            Epetra_MultiVector & C) const 
  {
    // do nothing, preconditioner is identity
    std::cout << "identiy precond" << std::endl;
  }
};

}} //namespace pressio::apps
#endif  // APPS_BURGERS1D_APPS_BURGERS1D_EPETRA_PRECONDITIONED_HPP_
