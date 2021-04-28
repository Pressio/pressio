/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_hypred_vel_problem_implicit_stepping.hpp
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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_HYPRED_VEL_PROBLEM_IMPLICIT_STEPPING_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_HYPRED_VEL_PROBLEM_IMPLICIT_STEPPING_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <typename ...Args>
class HypRedVeloProblemImplicitStepContinuousTimeApi
{

public:
  using this_t = HypRedVeloProblemImplicitStepContinuousTimeApi<Args...>;
  using traits = ::pressio::rom::details::traits<this_t>;

  using fom_system_t		= typename traits::fom_system_t;
  using scalar_t		= typename traits::scalar_t;
  using fom_native_state_t	= typename traits::fom_native_state_t;
  using fom_state_t		= typename traits::fom_state_t;
  using fom_velocity_t		= typename traits::fom_velocity_t;
  using galerkin_state_t	= typename traits::galerkin_state_t;
  using galerkin_native_state_t	= typename traits::galerkin_native_state_t;
  using galerkin_residual_t	= typename traits::galerkin_residual_t;
  using galerkin_jacobian_t	= typename traits::galerkin_jacobian_t;
  using decoder_t		= typename traits::decoder_t;
  using fom_state_reconstr_t	= typename traits::fom_state_reconstr_t;
  using fom_states_manager_t	= typename traits::fom_states_manager_t;
  using ud_ops_t		= typename traits::ud_ops_t;
  using projector_t		= typename traits::projector_t;
  using residual_policy_t	= typename traits::residual_policy_t;
  using jacobian_policy_t	= typename traits::jacobian_policy_t;
  using aux_stepper_t		= typename traits::aux_stepper_t;
  using stepper_t		= typename traits::stepper_t;
  static constexpr auto binding_sentinel = traits::binding_sentinel;

private:
  using At  = ::pressio::rom::impl::FomObjMixin<fom_system_t, binding_sentinel>;
  using Bt  = ::pressio::rom::impl::FomStatesMngrMixin<At, ud_ops_t, fom_state_t,
						       fom_state_reconstr_t, fom_states_manager_t>;
  using Ct  = HypRedVeloImplicitPoliciesMixin<Bt, residual_policy_t, jacobian_policy_t>;
  using m_t = ::pressio::rom::impl::ImplicitStepperMixin<Ct, aux_stepper_t, stepper_t>;
  m_t members_;

public:
  stepper_t & stepperRef(){ return members_.stepperObj_; }

  const fom_native_state_t & currentFomStateCRef() const{
    return *(members_.fomStatesMngr_.currentFomStateCRef().data());
  }

  const fom_state_reconstr_t & fomStateReconstructorCRef() const{
    return members_.fomStateReconstructor_;
  }

public:
  HypRedVeloProblemImplicitStepContinuousTimeApi() = delete;
  HypRedVeloProblemImplicitStepContinuousTimeApi(const HypRedVeloProblemImplicitStepContinuousTimeApi &) = default;
  HypRedVeloProblemImplicitStepContinuousTimeApi & operator=(const HypRedVeloProblemImplicitStepContinuousTimeApi &) = delete;
  HypRedVeloProblemImplicitStepContinuousTimeApi(HypRedVeloProblemImplicitStepContinuousTimeApi &&) = default;
  HypRedVeloProblemImplicitStepContinuousTimeApi & operator=(HypRedVeloProblemImplicitStepContinuousTimeApi &&) = delete;
  ~HypRedVeloProblemImplicitStepContinuousTimeApi() = default;

  /* ud_ops_t = void */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<std::is_void<_ud_ops_t>::value, int> = 0
    >
  HypRedVeloProblemImplicitStepContinuousTimeApi(const fom_system_t       & fomObj,
						 decoder_t	  & decoder,
						 const galerkin_state_t   & romStateIn,
						 const fom_native_state_t & fomNominalStateNative,
						 const projector_t & projector)
    : members_(romStateIn, fomObj, decoder, fomNominalStateNative, projector){}

  /* ud_ops_t != void */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<mpl::not_void<_ud_ops_t>::value, int> = 0
    >
  HypRedVeloProblemImplicitStepContinuousTimeApi(const fom_system_t     & fomObj,
						decoder_t	       & decoder,
						const galerkin_state_t   & romStateIn,
						const fom_native_state_t & fomNominalStateNative,
						const projector_t & projector,
						const _ud_ops_t & udOps)
    : members_(romStateIn, fomObj, decoder, fomNominalStateNative, projector, udOps){}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  template <
    bool _binding_sentinel = binding_sentinel,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t< _binding_sentinel and std::is_void<_ud_ops_t>::value,
      int > = 0
    >
  HypRedVeloProblemImplicitStepContinuousTimeApi(pybind11::object fomObjPython,
						 decoder_t & decoder,
						 const galerkin_native_state_t & romStateIn,
						 const fom_native_state_t fomNominalStateIn,
						 const projector_t & projector)
    : members_(galerkin_state_t(romStateIn), fomObjPython, decoder,
	       fomNominalStateIn, projector)
  {}
#endif
};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_HYPRED_VEL_PROBLEM_IMPLICIT_STEPPING_HPP_
