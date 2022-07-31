/*
//@HEADER
// ************************************************************************
//
// solvers_default.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_DEFAULT_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_DEFAULT_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

class DefaultUpdater
{
public:
  DefaultUpdater() = default;
  DefaultUpdater(DefaultUpdater const &) = default;
  DefaultUpdater & operator=(DefaultUpdater const &) = default;
  DefaultUpdater(DefaultUpdater &&) = default;
  DefaultUpdater & operator=(DefaultUpdater &&) = default;
  ~DefaultUpdater() = default;

public:
  void reset(){}

  template<typename SystemType, typename StateType, typename SolverMixinType>
  void operator()(const SystemType & /*unused*/,
		  StateType & state,
		  SolverMixinType & solver)
  {

    PRESSIOLOG_DEBUG("nonlinsolver: default update");
    // default update: y = y + alpha*correction
    const auto & correction = solver.correctionCRef();
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::update(state, one, correction, one);
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_DEFAULT_HPP_
