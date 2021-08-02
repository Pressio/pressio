/*
//@HEADER
// ************************************************************************
//
// solvers_create_gauss_newton.hpp
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

#ifndef SOLVERS_NONLINEAR_SOLVERS_CREATE_SOLVER_API_HPP_
#define SOLVERS_NONLINEAR_SOLVERS_CREATE_SOLVER_API_HPP_

#include "./impl/solvers_nonlinear_compose.hpp"

namespace pressio{ namespace nonlinearsolvers{

template<typename SystemType, typename StateType, typename ...Args>
auto createNewtonRaphson(const SystemType & system,
       const StateType & state,
       Args && ...args)
  -> impl::composeNewtonRaphson_t<SystemType, Args...>
{
  return impl::composeNewtonRaphson_t<SystemType, Args...>
    (system, state, std::forward<Args>(args)...);
}

//***** GN with NEQ *******
template<typename SystemType, typename StateType, typename ...Args>
auto createGaussNewton(const SystemType & system,
		       const StateType & state,
		       Args && ... args)
  -> impl::composeGaussNewton_t<SystemType, Args...>
{
  return impl::composeGaussNewton_t<SystemType, Args...>
    (system, state, std::forward<Args>(args)...);
}

//***** GN with QR *******
template<typename SystemType, typename StateType, typename ...Args>
auto createGaussNewtonQR(const SystemType & system,
			 const StateType & state,
			 Args && ...args)
  -> impl::composeGaussNewtonQR_t<SystemType, Args...>
{
  return impl::composeGaussNewtonQR_t<SystemType, Args...>
    (system, state, std::forward<Args>(args)...);
}

//***** IRWGN *******
namespace experimental{
template<typename SystemType, typename StateType, typename linear_solver_t>
auto createIRWGaussNewton(const SystemType & system,
			  const StateType & state,
			  linear_solver_t && linSolver)
  -> impl::composeIrwGaussNewton_t<SystemType, linear_solver_t>
{
  using c_t = impl::composeIrwGaussNewton<SystemType, linear_solver_t>;
  using w_t = typename c_t::weighting_t;
  using return_t = typename c_t::type;

  w_t W(system);
  return return_t(system, state,
		  std::forward<linear_solver_t>(linSolver),
		  std::move(W));
}
}// end namespace experimental

template<typename SystemType, typename StateType, typename ...Args>
auto createLevenbergMarquardt(const SystemType & system,
            const StateType & state,
            Args && ...args)
  -> impl::composeLevenbergMarquardt_t<SystemType, Args...>
{
  return impl::composeLevenbergMarquardt_t<SystemType, Args...>
  ( system, state, std::forward<Args>(args)...);
}

}}
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_GAUSS_NEWTON_HPP_