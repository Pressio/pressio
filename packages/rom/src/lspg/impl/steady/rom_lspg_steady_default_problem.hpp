/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_default_problem.hpp
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

#ifndef ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_DEFAULT_PROBLEM_HPP_
#define ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_DEFAULT_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace steady{

template <typename ...Args>
class DefaultProblemSteady
{
public:
  using this_t = DefaultProblemSteady<Args...>;
  using traits = ::pressio::rom::details::traits<this_t>;

  using fom_system_t		= typename traits::fom_system_t;
  using fom_native_state_t	= typename traits::fom_native_state_t;
  using fom_state_t		= typename traits::fom_state_t;
  using lspg_state_t		= typename traits::lspg_state_t;
  using decoder_t		= typename traits::decoder_t;
  using fom_state_reconstr_t	= typename traits::fom_state_reconstr_t;
  using fom_states_manager_t	= typename traits::fom_states_manager_t;
  using lspg_matrix_t		= typename traits::lspg_matrix_t;
  using residual_policy_t	= typename traits::residual_policy_t;
  using jacobian_policy_t	= typename traits::jacobian_policy_t;
  using system_t		= typename traits::system_t;

private:
  fom_state_t		fomStateReference_;
  fom_state_reconstr_t	fomStateReconstructor_;
  fom_states_manager_t	fomStatesMngr_;
  residual_policy_t	residualPolicy_;
  jacobian_policy_t	jacobianPolicy_;
  system_t		systemObj_;

public:
  system_t & systemRef(){
    return systemObj_;
  }

  const fom_state_reconstr_t & fomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

  const fom_native_state_t & currentFomStateCRef() const{
    return *fomStatesMngr_.currentFomStateCRef().data();
  }

public:
  DefaultProblemSteady() = delete;
  DefaultProblemSteady(const DefaultProblemSteady &) = default;
  DefaultProblemSteady & operator=(const DefaultProblemSteady &) = default;
  DefaultProblemSteady(DefaultProblemSteady &&) = default;
  DefaultProblemSteady & operator=(DefaultProblemSteady &&) = default;
  ~DefaultProblemSteady() = default;

  /* specialize for when the fom_system_t is regular c++ */
  template <
    typename _fom_system_t = fom_system_t,
    ::pressio::mpl::enable_if_t<
      !::pressio::ops::predicates::is_object_pybind<_fom_system_t>::value,
      int> = 0
  >
  DefaultProblemSteady(const _fom_system_t & fomSystemObj,
		       const decoder_t	& decoder,
		       const lspg_state_t & romStateIn,
		       const fom_native_state_t & fomNativeReferenceState)
    : fomStateReference_(fomNativeReferenceState),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      residualPolicy_(fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, decoder),
      systemObj_(fomSystemObj, residualPolicy_, jacobianPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /* specialize for when the fom_system_t is python object */
  template <
    typename _fom_system_t = fom_system_t,
    typename _lspg_state_t = lspg_state_t,
    ::pressio::mpl::enable_if_t<
      ::pressio::ops::predicates::is_object_pybind<_fom_system_t>::value and
      ::pressio::containers::predicates::is_vector_wrapper_pybind<_lspg_state_t>::value,
      int > = 0
  >
  DefaultProblemSteady(const _fom_system_t & fomSystemObj,
		       const decoder_t	& decoder,
		       typename ::pressio::containers::details::traits<_lspg_state_t>::wrapped_t & romStateIn,
		       const fom_native_state_t fomNativeReferenceState)
    : fomStateReference_(fomNativeReferenceState),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      residualPolicy_(fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, decoder),
      systemObj_(fomSystemObj, residualPolicy_, jacobianPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(_lspg_state_t(romStateIn));
  }
#endif

};

}}}}}
#endif  // ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_DEFAULT_PROBLEM_HPP_
