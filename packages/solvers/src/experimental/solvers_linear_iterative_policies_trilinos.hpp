/*
//@HEADER
// ************************************************************************
//
// solvers_linear_iterative_policies_trilinos.hpp
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

#ifdef HAVE_TRILINOS
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_POLICIES_TRILINOS_HPP_
#define SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_POLICIES_TRILINOS_HPP_

#include "AztecOO.h"

namespace pressio{
namespace solvers{
namespace trilinos_policies {


// Linear solvers policies
struct CG {
  static void set_solver_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_solver, AZ_cg);
  }
};


struct Gmres {
  static void set_solver_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_solver, AZ_gmres);
  }
};


struct BicgStab {
  static void set_solver_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_solver, AZ_bicgstab);
  }
};


// Preconditioners policies
struct DefaultPreconditioner {
  static void set_preconditioner_type(AztecOO& solver) {
    return;
  }
};

struct DecompIlu {
  static void set_preconditioner_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
  }
};


struct DecompIlut {
  static void set_preconditioner_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    solver.SetAztecOption(AZ_overlap, 1);
    solver.SetAztecOption(AZ_ilut_fill, 3.0);
  }
};


struct DecompIcc {
  static void set_preconditioner_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve, AZ_icc); 
  }
};


struct Jacobi {
  static void set_preconditioner_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_precond, AZ_Jacobi);
    solver.SetAztecOption(AZ_poly_ord, 1);
  }
};


} // end namespace trilinos_policies
} // end namespace solvers
}//end namespace pressio
#endif
#endif
