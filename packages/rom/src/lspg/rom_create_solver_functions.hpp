/*
//@HEADER
// ************************************************************************
//
// rom_create_solver_functions.hpp
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

#ifndef ROM_LSPG_ROM_CREATE_SOLVER_FUNCTIONS_HPP_
#define ROM_LSPG_ROM_CREATE_SOLVER_FUNCTIONS_HPP_

#include "./impl/rom_view_problem_system.hpp"

namespace pressio{ namespace rom{ namespace lspg{

/* why using this over just instantiating a solver directly?
   (1) we can limit the solver types admissible for a galerkin problem
   we know that right now only newton-raphson makes sense for galerkin

   (2) allows users to pass a rom problem object without needing
   to know that the rom problem has a stepper object inside which would
   need to be extracted and passed to the solver

   (3) potentially, we could specialize these functions to pick the
   best solver/parameters ourselves based on some conditions
*/

/* ========================
    Gauss-Newton normal-eq
   ========================*/
template<typename rom_problem_t, typename ...Args>
auto createGaussNewtonSolver(rom_problem_t & problem, Args && ... args)
->
  decltype(::pressio::solvers::nonlinear::createGaussNewton
    (impl::_SystemOrStepper(problem), std::forward<Args>(args)...))
{
  return ::pressio::solvers::nonlinear::createGaussNewton
    (impl::_SystemOrStepper(problem), std::forward<Args>(args)...);
}

/* ========================
	Gauss-Newton QR
   ========================*/
template<typename rom_problem_t, typename ...Args>
auto createGaussNewtonQRSolver(rom_problem_t & problem, Args && ... args)
->  decltype(::pressio::solvers::nonlinear::createGaussNewtonQR
    (impl::_SystemOrStepper(problem), std::forward<Args>(args)...))
{
  return ::pressio::solvers::nonlinear::createGaussNewtonQR
    (impl::_SystemOrStepper(problem), std::forward<Args>(args)...);
}

/* ========================
      LevenbergMarquardt
   ========================*/
template<typename rom_problem_t, typename ...Args>
auto createLevenbergMarquardtSolver(rom_problem_t & problem, Args && ... args)
->  decltype(::pressio::solvers::nonlinear::createLevenbergMarquardt
    (impl::_SystemOrStepper(problem), std::forward<Args>(args)...))
{
  return ::pressio::solvers::nonlinear::createLevenbergMarquardt
    (impl::_SystemOrStepper(problem), std::forward<Args>(args)...);
}

}}}//end namespace pressio::rom::lspg
#endif  // ROM_LSPG_ROM_CREATE_SOLVER_FUNCTIONS_HPP_
