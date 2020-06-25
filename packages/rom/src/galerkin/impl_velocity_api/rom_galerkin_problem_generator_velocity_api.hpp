/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_problem_generator.hpp
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

#ifndef ROM_GALERKIN_PROBLEM_GENERATOR_VELOCITY_API_HPP_
#define ROM_GALERKIN_PROBLEM_GENERATOR_VELOCITY_API_HPP_


namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <
  template <class ...> class galerkin_type,
  typename fom_type,
  typename stepper_tag,
  typename rom_state_type,
  typename ...Args
  >
class ProblemGeneratorVelocityApi
{

public:
  // define the type holding types for the problem
  using problem_t = galerkin_type<fom_type, stepper_tag, rom_state_type, Args...>;

  using fom_t			= typename problem_t::fom_t;
  using scalar_t		= typename problem_t::scalar_t;
  using fom_native_state_t	= typename problem_t::fom_native_state_t;
  using fom_state_t		= typename problem_t::fom_state_t;
  using fom_velocity_t		= typename problem_t::fom_velocity_t;

  using galerkin_state_t	= typename problem_t::galerkin_state_t;
  using galerkin_native_state_t	= typename problem_t::galerkin_native_state_t;
  using decoder_t		= typename problem_t::decoder_t;
  using fom_state_reconstr_t	= typename problem_t::fom_state_reconstr_t;
  using fom_states_manager_t		= typename problem_t::fom_states_manager_t;
  using ud_ops_t		= typename problem_t::ud_ops_t;

  using residual_policy_t	= typename problem_t::residual_policy_t;
  using stepper_t		= typename problem_t::stepper_t;

private:
  fom_state_t			fomStateReference_;
  fom_state_reconstr_t		fomStateReconstructor_;
  fom_velocity_t		fomVelocityRef_;
  fom_states_manager_t		fomStatesMngr_;
  residual_policy_t		residualPolicy_;
  stepper_t			stepperObj_;

public:
  stepper_t & getStepperRef(){
    return stepperObj_;
  }

  const fom_state_reconstr_t & getFomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

public:
  ProblemGeneratorVelocityApi() = delete;
  ~ProblemGeneratorVelocityApi() = default;

  /*
   * ud_ops_t = void, C++ types
  */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops_t>::value and
      ::pressio::containers::meta::is_wrapper<galerkin_state_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !::pressio::containers::meta::is_vector_wrapper_pybind<galerkin_state_t>::value
#endif
      , int> = 0
  >
  ProblemGeneratorVelocityApi(const fom_t   & appObj,
			      const fom_native_state_t & yFomRefNative,
			      const decoder_t	    & decoder,
			      galerkin_state_t	    & yROM,
			      scalar_t		    t0)
    : fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomVelocityRef_( appObj.velocity(*fomStateReference_.data(), t0) ),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      residualPolicy_(fomVelocityRef_, fomStatesMngr_, decoder),
      stepperObj_(yROM, appObj, residualPolicy_)
  {}

  /*
   * ud_ops_t != void, C++ types
  */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      !std::is_void<_ud_ops_t>::value and
      ::pressio::containers::meta::is_wrapper<galerkin_state_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !::pressio::containers::meta::is_vector_wrapper_pybind<galerkin_state_t>::value
#endif
      , int> = 0
  >
  ProblemGeneratorVelocityApi(const fom_t   & appObj,
			      const fom_native_state_t & yFomRefNative,
			      const decoder_t	    & decoder,
			      galerkin_state_t	    & yROM,
			      scalar_t		    t0,
			      const _ud_ops_t & udOps)
    : fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder, udOps),
      fomVelocityRef_( appObj.velocity(*fomStateReference_.data(), t0) ),
      fomStatesMngr_(fomStateReconstructor_, &udOps, fomStateReference_),
      residualPolicy_(fomVelocityRef_, fomStatesMngr_, decoder, udOps),
      stepperObj_(yROM, appObj, residualPolicy_)
  {}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /*
   * ud_ops_t == void and state_type is wrapper of pybind11::array
  */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops_t>::value and
      ::pressio::containers::meta::is_vector_wrapper_pybind<galerkin_state_t>::value, 
      int > = 0
  >
  ProblemGeneratorVelocityApi(const fom_t   & appObj,
			      fom_native_state_t	    yFomRefNative,
			      const decoder_t	    & decoder,
			      galerkin_native_state_t  yROM,
			      scalar_t		    t0)
    : fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomVelocityRef_( appObj.attr("velocity")(*fomStateReference_.data(), t0) ),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      residualPolicy_(fomVelocityRef_, fomStatesMngr_, decoder),
      stepperObj_(galerkin_state_t(yROM), appObj, residualPolicy_)
  {}
#endif
};

}}}}//end namespace pressio::rom::galerkin::impl
#endif
