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

private:
  fom_state_t			fomStateReference_;
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
  DefaultProblemContinuousTimeApi & operator=(const DefaultProblemContinuousTimeApi &) = default;
  DefaultProblemContinuousTimeApi(DefaultProblemContinuousTimeApi &&) = default;
  DefaultProblemContinuousTimeApi & operator=(DefaultProblemContinuousTimeApi &&) = default;
  ~DefaultProblemContinuousTimeApi() = default;

  /*
   * ud_ops_t = void, C++ types
   */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops_t>::value and
      ::pressio::containers::predicates::is_wrapper<galerkin_state_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !::pressio::containers::predicates::is_vector_wrapper_pybind<galerkin_state_t>::value
#endif
      , int> = 0
    >
  DefaultProblemContinuousTimeApi(const fom_system_t	   & fomSystemObj,
				  const decoder_t	   & decoder,
				  const galerkin_state_t   & romStateIn,
				  const fom_native_state_t & fomNativeReferenceState)
    : fomStateReference_(fomNativeReferenceState),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomVelocityRef_(fomSystemObj.createVelocity()),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      velocityPolicy_(fomVelocityRef_, fomStatesMngr_, decoder),
      stepperObj_(romStateIn, fomSystemObj, velocityPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

  /*
   * ud_ops_t != void, C++ types
   */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      !std::is_void<_ud_ops_t>::value and
      ::pressio::containers::predicates::is_wrapper<galerkin_state_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !::pressio::containers::predicates::is_vector_wrapper_pybind<galerkin_state_t>::value
#endif
      , int> = 0
    >
  DefaultProblemContinuousTimeApi(const fom_system_t	    & fomSystemObj,
				  const decoder_t	    & decoder,
				  const galerkin_state_t    & romStateIn,
				  const fom_native_state_t  & fomNativeReferenceState,
				  const _ud_ops_t	    & udOps)
    : fomStateReference_(fomNativeReferenceState),
      fomStateReconstructor_(fomStateReference_, decoder, udOps),
      fomVelocityRef_(fomSystemObj.createVelocity()),
      fomStatesMngr_(fomStateReconstructor_, &udOps, fomStateReference_),
      velocityPolicy_(fomVelocityRef_, fomStatesMngr_, decoder, udOps),
      stepperObj_(romStateIn, fomSystemObj, velocityPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /*
   * ud_ops_t == void and state_type is wrapper of pybind11::array
   */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops_t>::value and
      ::pressio::containers::predicates::is_vector_wrapper_pybind<galerkin_state_t>::value,
      int > = 0
    >
  DefaultProblemContinuousTimeApi(const fom_system_t	    & fomSystemObj,
				  const fom_native_state_t & fomNativeReferenceState,
				  const decoder_t	    & decoder,
				  const galerkin_native_state_t & romStateIn,
				  scalar_t		    t0)
    : fomStateReference_(fomNativeReferenceState),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomVelocityRef_( fomSystemObj.attr("velocity")(*fomStateReference_.data(), t0) ),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      velocityPolicy_(fomVelocityRef_, fomStatesMngr_, decoder),
      stepperObj_(galerkin_state_t(romStateIn), fomSystemObj, velocityPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(galerkin_state_t(romStateIn));
  }
#endif

};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_CONTINUOUS_TIME_API_HPP_
