/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_default_problem_continuous_time_api.hpp
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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_CONTINUOUS_TIME_API_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_CONTINUOUS_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template<bool forPy, typename T> struct _systemMemberType;

template<typename T> struct _systemMemberType<true, T> {
  using type = T; };
template<typename T> struct _systemMemberType<false, T>{
  using type = ::pressio::utils::impl::empty;};


template <typename ...Args>
class DefaultProblemContinuousTimeApi
{

public:
  using this_t = DefaultProblemContinuousTimeApi<Args...>;
  using traits = ::pressio::rom::details::traits<this_t>;

  using fom_system_t		= typename traits::fom_system_t;
  using scalar_t		= typename traits::scalar_t;
  using fom_native_state_t	= typename traits::fom_native_state_t;
  using fom_state_t		= typename traits::fom_state_t;
  using fom_velocity_t		= typename traits::fom_velocity_t;

  using galerkin_state_t	= typename traits::galerkin_state_t;
  using galerkin_native_state_t	= typename traits::galerkin_native_state_t;
  using decoder_t		= typename traits::decoder_t;
  using fom_state_reconstr_t	= typename traits::fom_state_reconstr_t;
  using fom_states_manager_t	= typename traits::fom_states_manager_t;
  using ud_ops_t		= typename traits::ud_ops_t;

  using velocity_policy_t	= typename traits::velocity_policy_t;
  using stepper_t		= typename traits::stepper_t;

  static constexpr auto binding_sentinel = traits::binding_sentinel;

private:
  // when dealing with pressio4py, the fom_system_t is a C++ class in pressio4py
  // that wraps the actual FOM python object. to construct this ROM problem,
  // the Python code passes the python FOM object NOT a C++ object instantiated
  // from fom_system_t, see constructor at the end for pressio4py.
  // Therefore, ONLY when we deal with pressio4py, we make this problem store
  // an object of the the fom_system_t.
  typename _systemMemberType<binding_sentinel, fom_system_t>::type fomObj_;

  fom_state_t			fomNominalState_;
  fom_state_reconstr_t		fomStateReconstructor_;
  fom_velocity_t		fomVelocityRef_;
  fom_states_manager_t		fomStatesMngr_;
  velocity_policy_t		velocityPolicy_;
  stepper_t			stepperObj_;

public:
  stepper_t & stepperRef(){
    return stepperObj_;
  }

  const fom_state_reconstr_t & fomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

  const fom_native_state_t & currentFomStateCRef() const{
    return *fomStatesMngr_.currentFomStateCRef().data();
  }

public:
  DefaultProblemContinuousTimeApi() = delete;
  DefaultProblemContinuousTimeApi(const DefaultProblemContinuousTimeApi &) = default;
  DefaultProblemContinuousTimeApi & operator=(const DefaultProblemContinuousTimeApi &) = delete;
  DefaultProblemContinuousTimeApi(DefaultProblemContinuousTimeApi &&) = default;
  DefaultProblemContinuousTimeApi & operator=(DefaultProblemContinuousTimeApi &&) = delete;
  ~DefaultProblemContinuousTimeApi() = default;

  /*
   * ud_ops_t = void, not binding
   */
  template <
    bool _binding_sentinel=binding_sentinel,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops_t>::value and
      _binding_sentinel==false, int> = 0
    >
  DefaultProblemContinuousTimeApi(const fom_system_t	   & fomObj,
				  const decoder_t	   & decoder,
				  const galerkin_state_t   & romStateIn,
				  const fom_native_state_t & fomNominalStateNative)
    : fomNominalState_(fomNominalStateNative),
      fomStateReconstructor_(fomNominalState_, decoder),
      fomVelocityRef_(fomObj.createVelocity()),
      fomStatesMngr_(fomStateReconstructor_, fomNominalState_),
      velocityPolicy_(fomVelocityRef_, fomStatesMngr_, decoder),
      stepperObj_(romStateIn, fomObj, velocityPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

  /*
   * ud_ops_t != void, not binding
   */
  template <
    bool _binding_sentinel=binding_sentinel,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      !std::is_void<_ud_ops_t>::value and
      _binding_sentinel==false, int> = 0
    >
  DefaultProblemContinuousTimeApi(const fom_system_t	    & fomObj,
				  const decoder_t	    & decoder,
				  const galerkin_state_t    & romStateIn,
				  const fom_native_state_t  & fomNominalStateNative,
				  const _ud_ops_t	    & udOps)
    : fomNominalState_(fomNominalStateNative),
      fomStateReconstructor_(fomNominalState_, decoder, udOps),
      fomVelocityRef_(fomObj.createVelocity()),
      fomStatesMngr_(fomStateReconstructor_, fomNominalState_, udOps),
      velocityPolicy_(fomVelocityRef_, fomStatesMngr_, decoder, udOps),
      stepperObj_(romStateIn, fomObj, velocityPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /* ud_ops_t == void and state_type is wrapper of pybind11::array

     Note that in this case, the constructor receives a pybind11::object
     as the FOM object, and we use that to construct the C++ wrapper object.
  */
  template <
    bool _binding_sentinel=binding_sentinel,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops_t>::value and
      _binding_sentinel, int> = 0
    >
  DefaultProblemContinuousTimeApi(pybind11::object fomObjPython,
				  const decoder_t & decoder,
				  const galerkin_native_state_t & romStateIn,
				  const fom_native_state_t & fomNominalStateNative)
    : fomObj_(fomObjPython),
      fomNominalState_(fomNominalStateNative),
      fomStateReconstructor_(fomNominalState_, decoder),
      fomVelocityRef_(fomObj_.createVelocity()),
      fomStatesMngr_(fomStateReconstructor_, fomNominalState_),
      velocityPolicy_(fomVelocityRef_, fomStatesMngr_, decoder),
      stepperObj_(galerkin_state_t(romStateIn), fomObj_, velocityPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    auto romStateView = galerkin_state_t(romStateIn, ::pressio::view());
    fomStatesMngr_.reconstructCurrentFomState(romStateView);
  }
#endif
};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_CONTINUOUS_TIME_API_HPP_