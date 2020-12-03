/*
//@HEADER
// ************************************************************************
//
// solvers_jacobian_action_wrapper.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_JACOBIAN_ACTION_WRAPPER_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_JACOBIAN_ACTION_WRAPPER_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<class state_type, class system_type>
struct JacobianActionWrapper
{
  state_type const & state_;
  system_type const & sys_;
  bool updateJacobianAction_ = true;

  // for now, delete all specials mf since we need to fix first the
  // memebers to not be just simple references but someting else
  // that have proper semantics
  JacobianActionWrapper() = delete;
  JacobianActionWrapper(JacobianActionWrapper const &) = delete;
  JacobianActionWrapper & operator=(JacobianActionWrapper const &) = delete;
  JacobianActionWrapper(JacobianActionWrapper && o) = delete;
  JacobianActionWrapper & operator=(JacobianActionWrapper && o) = delete;
  ~JacobianActionWrapper() = default;

  JacobianActionWrapper(state_type const & state,
	      system_type const & system,
	      bool updateJacobianAction)
    : state_(state), sys_(system),
      updateJacobianAction_(updateJacobianAction)
  {}

  template<typename T>
  void applyJacobian(const T & in, T & out) const
  {
    sys_.applyJacobian(state_, in, out, updateJacobianAction_);
  }
};

}//end namespace impl

template<typename ... Args>
using JacobianActionWrapper = impl::JacobianActionWrapper<Args...>;

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_JACOBIAN_ACTION_WRAPPER_HPP_
