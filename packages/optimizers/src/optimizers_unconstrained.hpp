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

#ifndef OPTIMIZERS_UNCONSTRAINED_HPP_
#define OPTIMIZERS_UNCONSTRAINED_HPP_

#include "ROL_OptimizationSolver.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StdObjective.hpp"

namespace pressio{ namespace optimizers{

namespace impl{

template<typename scalar_type, typename system_type>
class RolObjectiveWrapper : public ROL::StdObjective<scalar_type>
{
  using state_type = typename system_type::state_type;
  const system_type & wrappedObj_;

public:
  RolObjectiveWrapper(const system_type & wrappedObj) : wrappedObj_(wrappedObj){}

  scalar_type value(const std::vector<scalar_type> &x, scalar_type &tol){
    state_type x2(x.size());
    for (auto i=0; i<x2.extent(0); ++i) x2[i] = x[i];

    return wrappedObj_(x2);
  }
};

template <typename system_type>
class UnconstrainedRol
{
  using scalar_t = typename system_type::scalar_type;
  using state_t  = typename system_type::state_type;

  const ::pressio::optimizers::Parameters<scalar_t> & params_;
  ROL::ParameterList rolParList_;

public:
  UnconstrainedRol(const ::pressio::optimizers::Parameters<scalar_t> & params) : params_(params)
  {
    // convert params to rolParList
    rolParList_.sublist("Step").set("Type","Trust Region");
    rolParList_.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
  }

  void solve(const system_type & sysObj, state_t & optState)
  {
    using wrapper_t = RolObjectiveWrapper<scalar_t, system_type>;
    auto obj = ROL::makePtr<wrapper_t>(sysObj);

    const auto optSize = optState.extent(0);
    ROL::Ptr<ROL::StdVector<scalar_t> > x = ROL::makePtr<ROL::StdVector<scalar_t>>(optSize);
    for (auto i=0; i<optState.extent(0); ++i) optState[i] = (*x)[i];

    ROL::OptimizationProblem<scalar_t> problem( obj, x );
    problem.check(std::cout);
    ROL::OptimizationSolver<scalar_t> solver( problem, rolParList_ );
    solver.solve(std::cout);
  }
};

} //end namespace impl

template <typename ...Args>
using Unconstrained = impl::UnconstrainedRol<Args...>;

}}//end namespace pressio::optimizers
#endif
