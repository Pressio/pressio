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

template<bool forPy, typename T> struct _systemMemberType;

template<typename T> struct _systemMemberType<true, T> {
  using type = T; };
template<typename T> struct _systemMemberType<false, T>{
  using type = ::pressio::utils::impl::empty;};


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
  static constexpr auto binding_sentinel = traits::binding_sentinel;
  using lspg_native_state_t =
    typename ::pressio::containers::details::traits<lspg_state_t>::wrapped_t;

private:
  // when dealing with pressio4py, the fom_system_t is a C++ class in pressio4py
  // that wraps the actual FOM python object. to construct this ROM problem,
  // the Python code passes the python FOM object NOT a C++ object instantiated
  // from fom_system_t, see constructor at the end for pressio4py.
  // Therefore, ONLY when we deal with pressio4py, we make this problem store
  // an object of the fom_system_t.
  typename _systemMemberType<binding_sentinel, fom_system_t>::type fomObj_;

  fom_state_t		fomNominalState_;
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

  template <
    bool _binding_sentinel = binding_sentinel,
    ::pressio::mpl::enable_if_t<!_binding_sentinel, int> = 0
    >
  DefaultProblemSteady(const fom_system_t & fomObj,
		       const decoder_t	& decoder,
		       const lspg_state_t & romStateIn,
		       const fom_native_state_t & fomNominalNative)
    : fomNominalState_(fomNominalNative),
      fomStateReconstructor_(fomNominalState_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomNominalState_),
      residualPolicy_(fomStatesMngr_),
      jacobianPolicy_(fomStatesMngr_, decoder),
      systemObj_(fomObj, residualPolicy_, jacobianPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  template <
    bool _binding_sentinel = binding_sentinel,
    ::pressio::mpl::enable_if_t<_binding_sentinel, int> = 0
    >
  DefaultProblemSteady(pybind11::object fomObjPython,
		       const decoder_t & decoder,
		       const lspg_native_state_t & romStateIn,
		       const fom_native_state_t fomNominalNative)
  : fomObj_(fomObjPython),
    fomNominalState_(fomNominalNative),
    fomStateReconstructor_(fomNominalState_, decoder),
    fomStatesMngr_(fomStateReconstructor_, fomNominalState_),
    residualPolicy_(fomStatesMngr_),
    jacobianPolicy_(fomStatesMngr_, decoder),
    systemObj_(fomObj_, residualPolicy_, jacobianPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(lspg_state_t(romStateIn));
  }
#endif
};

}}}}}
#endif  // ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_DEFAULT_PROBLEM_HPP_
