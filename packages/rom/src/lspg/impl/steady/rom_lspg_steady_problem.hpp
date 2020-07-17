/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_problem.hpp
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

#ifndef ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_PROBLEM_HPP_
#define ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_PROBLEM_HPP_

#include "./policies/rom_lspg_steady_residual_policy.hpp"
#include "./policies/rom_lspg_steady_jacobian_policy.hpp"
#include "rom_lspg_steady_system.hpp"

#include "./traits/rom_lspg_steady_common_traits.hpp"
#include "./traits/rom_lspg_steady_default_problem_traits.hpp"
#include "./traits/rom_lspg_steady_preconditioned_problem_traits.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace steady{

template <
  template <class ...> class lspg_type,
  typename system_type,
  typename lspg_state_type,
  typename decoder_type,
  typename ...Args
  >
class ProblemSteady
{

public:
  // define the type holding types for the problem
  using lspg_problem_t = lspg_type<system_type, lspg_state_type, decoder_type, Args...>;

  using fom_t = typename lspg_problem_t::fom_t;
  using fom_native_state_t = typename lspg_problem_t::fom_native_state_t;
  using fom_state_t = typename lspg_problem_t::fom_state_t;
  using lspg_state_t = typename lspg_problem_t::lspg_state_t;
  using decoder_t = typename lspg_problem_t::decoder_t;
  using fom_state_reconstr_t = typename lspg_problem_t::fom_state_reconstr_t;
  using fom_states_manager_t = typename lspg_problem_t::fom_states_manager_t;
  using lspg_matrix_t = typename lspg_problem_t::lspg_matrix_t;
  using lspg_residual_policy_t = typename lspg_problem_t::lspg_residual_policy_t;
  using lspg_jacobian_policy_t = typename lspg_problem_t::lspg_jacobian_policy_t;
  using lspg_system_t = typename lspg_problem_t::lspg_system_t;

private:
  fom_state_t			fomStateReference_;
  fom_state_reconstr_t		fomStateReconstructor_;
  fom_states_manager_t		fomStatesMngr_;
  lspg_residual_policy_t	residualPolicy_;
  lspg_jacobian_policy_t	jacobianPolicy_;
  lspg_system_t			systemObj_;

public:
  lspg_system_t & getSystemRef(){
    return systemObj_;
  }

  const fom_state_reconstr_t & getFomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

public:
  /* specialize for when the fom_t is regular c++ */
  template <
    typename _fom_t = fom_t,
    ::pressio::mpl::enable_if_t<
      !::pressio::ops::predicates::is_object_pybind<_fom_t>::value,
      int> = 0
  >
  ProblemSteady(const _fom_t & appObj,
  		const fom_native_state_t & yFomRefNative,
  		const decoder_t	& decoder)
    : fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      residualPolicy_(fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, decoder),
      systemObj_(appObj, residualPolicy_, jacobianPolicy_)
  {}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /* specialize for when the fom_t is python object */
  template <
    typename _fom_t = fom_t,
    ::pressio::mpl::enable_if_t<
      ::pressio::ops::predicates::is_object_pybind<_fom_t>::value,
      int > = 0
  >
  ProblemSteady(const _fom_t & appObj,
		const fom_native_state_t yFomRefNative,
		const decoder_t	& decoder)
    : fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      residualPolicy_(fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, decoder),
      systemObj_(appObj, residualPolicy_, jacobianPolicy_)
  {}
#endif

};

}}}}}
#endif  // ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_PROBLEM_HPP_
