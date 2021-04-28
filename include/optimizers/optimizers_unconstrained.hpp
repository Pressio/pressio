/*
//@HEADER
// ************************************************************************
//
// optimizers_unconstrained.hpp
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

#ifndef OPTIMIZERS_OPTIMIZERS_UNCONSTRAINED_HPP_
#define OPTIMIZERS_OPTIMIZERS_UNCONSTRAINED_HPP_

#include "ROL_OptimizationSolver.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StdObjective.hpp"

namespace pressio{ namespace optimizers{

namespace impl{

template<typename scalar_type, typename system_type>
class RolObjectiveWrapper : public ROL::Objective<scalar_type>
{
  using state_type = typename system_type::state_type;
  const system_type & wrappedObj_;

public:
  RolObjectiveWrapper(const system_type & wrappedObj) : wrappedObj_(wrappedObj){}

  scalar_type value(const ROL::Vector<scalar_type> &x, scalar_type &tol)
  {
    const state_type & ex = static_cast<const state_type&>(x);
    return wrappedObj_(ex);
  }

  void gradient( ROL::Vector<scalar_type> &g,
		 const ROL::Vector<scalar_type> &x,
		 scalar_type &tol )
  {
    const state_type & x2 = static_cast<const state_type&>(x);
    state_type	     & g2 = static_cast<state_type&>(g);
    wrappedObj_.gradient(x2, g2);
  }
};


template <typename system_type>
class UnconstrainedRol
{
  using scalar_t = typename system_type::scalar_type;
  using state_t  = typename system_type::state_type;

  ROL::ParameterList rolParList_;
  const ::pressio::optimizers::Parameters<scalar_t> & params_;


public:
  UnconstrainedRol(const ::pressio::optimizers::Parameters<scalar_t> & params)
    : params_(params)
  {
    ::pressio::optimizers::convertToRolParameterList(params, rolParList_);
    rolParList_.print(std::cout);
  }

  void solve(const system_type & sysObj, state_t & optState)
  {
    using wrapper_t = RolObjectiveWrapper<scalar_t, system_type>;
    auto obj = ROL::makePtr<wrapper_t>(sysObj);

    // make ptr from reference, does not make any new allocation
    auto x = ROL::makePtrFromRef(optState);

    ROL::OptimizationProblem<scalar_t> problem( obj, x);
    problem.check(std::cout);
    ROL::OptimizationSolver<scalar_t> solver( problem, rolParList_ );
    solver.solve(std::cout);
  }
};

} //end namespace impl


template <typename ...Args>
using Unconstrained = impl::UnconstrainedRol<Args...>;


}}//end namespace pressio::optimizers
#endif  // OPTIMIZERS_OPTIMIZERS_UNCONSTRAINED_HPP_
