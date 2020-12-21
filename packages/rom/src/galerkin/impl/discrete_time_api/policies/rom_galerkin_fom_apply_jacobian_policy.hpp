/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_jacobian_policy_discrete_time_api.hpp
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

#ifndef ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_APPLY_JACOBIAN_POLICY_HPP_
#define ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_APPLY_JACOBIAN_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <class fom_states_manager_t, class fom_apply_jac_type, class decoder_type>
class FomApplyJacobianPolicyDiscreteTimeApi
{
public:
  using data_type = fom_apply_jac_type;

private:
  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  std::reference_wrapper<const typename decoder_type::jacobian_type> phi_;
  mutable fom_apply_jac_type fomApplyJac_;

public:
  FomApplyJacobianPolicyDiscreteTimeApi() = delete;
  FomApplyJacobianPolicyDiscreteTimeApi(const FomApplyJacobianPolicyDiscreteTimeApi &) = default;
  FomApplyJacobianPolicyDiscreteTimeApi & operator=(const FomApplyJacobianPolicyDiscreteTimeApi &) = delete;
  FomApplyJacobianPolicyDiscreteTimeApi(FomApplyJacobianPolicyDiscreteTimeApi &&) = default;
  FomApplyJacobianPolicyDiscreteTimeApi & operator=(FomApplyJacobianPolicyDiscreteTimeApi &&) = delete;
  ~FomApplyJacobianPolicyDiscreteTimeApi() = default;

  template<typename fom_system_t>
  FomApplyJacobianPolicyDiscreteTimeApi(const fom_system_t & fomSystemObj,
					fom_states_manager_t & fomStatesMngr,
					const decoder_type & decoder)
    : fomStatesMngr_(fomStatesMngr),
      phi_(decoder.jacobianCRef()),
      fomApplyJac_(fomSystemObj.createApplyDiscreteTimeJacobianResult(*phi_.get().data()))
  {}

public:
  const fom_apply_jac_type & get() const{
    return fomApplyJac_;
  }

  // we have here n = 1 prev rom states
  template<
    typename galerkin_state_t,
    typename galerkin_prev_states_t,
    typename fom_system_t,
    typename scalar_t
  >
  mpl::enable_if_t< galerkin_prev_states_t::size()==1 >
  compute(const galerkin_state_t	& galerkinState,
	  const fom_system_t	& fomSystemObj,
	  const scalar_t		& time,
	  const galerkin_prev_states_t & galerkinPrevStates,
	  const scalar_t		& dt,
	  const ::pressio::ode::types::step_t & step) const
  {
    // here we assume that the residual policy already:
    // - reconstucted the previous states

    const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
    const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
    fomSystemObj.applyDiscreteTimeJacobian(step, time, dt,
					   *(phi_.get().data()),
					   *fomApplyJac_.data(),
					   *yn.data(),
					   *ynm1.data());
  }

  // we have here n = 2 prev rom states
  template<
    typename galerkin_state_t,
    typename galerkin_prev_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  mpl::enable_if_t< galerkin_prev_states_t::size()==2 >
  compute(const galerkin_state_t	& galerkinState,
	  const fom_system_t       & fomSystemObj,
	  const scalar_t		& time,
	  const galerkin_prev_states_t & galerkinPrevStates,
	  const scalar_t		& dt,
	  const ::pressio::ode::types::step_t & step) const
  {
    // here we assume that the residual policy already:
    // - reconstucted the previous states

    const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
    const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
    const auto & ynm2 = fomStatesMngr_.get().fomStatePrevStepCRef();
    fomSystemObj.applyDiscreteTimeJacobian(step, time, dt,
					   *(phi_.get().data()),
					   *fomApplyJac_.data(),
					   *yn.data(),
					   *ynm1.data());
  }
};


// template<
//   typename galerkin_jacobian_type,
//   typename fom_states_manager_t,
//   typename fom_apply_jacobian_ret_type,
//   typename decoder_type,
//   typename projector_t
//   >
// class JacobianPolicyDiscreteTimeApi
// {
// public:
//   using time_step_type = ::pressio::ode::types::step_t;

// public:
//   JacobianPolicyDiscreteTimeApi() = delete;
//   JacobianPolicyDiscreteTimeApi(const JacobianPolicyDiscreteTimeApi &) = default;
//   JacobianPolicyDiscreteTimeApi & operator=(const JacobianPolicyDiscreteTimeApi &) = delete;
//   JacobianPolicyDiscreteTimeApi(JacobianPolicyDiscreteTimeApi &&) = default;
//   JacobianPolicyDiscreteTimeApi & operator=(JacobianPolicyDiscreteTimeApi &&) = delete;
//   ~JacobianPolicyDiscreteTimeApi() = default;

//   template<typename fom_system_t>
//   JacobianPolicyDiscreteTimeApi(std::size_t romSize,
// 				const fom_system_t & fomSystemObj,
// 				fom_states_manager_t & fomStatesMngr,
// 				const projector_t & projector,
// 				const decoder_type & decoder)
//     : romSize_(romSize),
//       fomStatesMngr_(fomStatesMngr),
//       decoderObj_(decoder),
//       phi_(decoder.jacobianCRef()),
//       fomApplyJac_(fomSystemObj.createApplyDiscreteTimeJacobianResult(*phi_.get().data())),
//       projector_(projector)
//   {}

// public:
//   template <typename fom_system_t>
//   galerkin_jacobian_type create(const fom_system_t & fomSystemObj) const
//   {
//     return galerkin_jacobian_type(romSize_, romSize_);
//   }

//   template <
//     typename stepper_tag,
//     typename galerkin_state_t,
//     typename galerkin_prev_states_t,
//     typename fom_system_t,
//     typename scalar_t
//     >
//   void compute(const galerkin_state_t	& galerkinState,
// 	       const galerkin_prev_states_t	& galerkinPrevStates,
// 	       const fom_system_t	& fomSystemObj,
// 	       const scalar_t		& time,
// 	       const scalar_t		& dt,
// 	       const time_step_type	& step,
// 	       galerkin_jacobian_type	& galerkinJacobian) const
//   {
//     this->compute_impl(galerkinState, galerkinPrevStates, fomSystemObj,
// 		       time, dt, step, galerkinJacobian);
//   }

// private:
//   // we have here n = 1 prev rom states
//   template<
//     typename galerkin_state_t,
//     typename galerkin_prev_states_t,
//     typename fom_system_t,
//     typename scalar_t
//   >
//   mpl::enable_if_t< galerkin_prev_states_t::size()==1 >
//   compute_impl(const galerkin_state_t	& galerkinState,
// 	       const galerkin_prev_states_t & galerkinPrevStates,
// 	       const fom_system_t	& fomSystemObj,
// 	       const scalar_t		& time,
// 	       const scalar_t		& dt,
// 	       const time_step_type	& step,
// 	       galerkin_jacobian_type	& galerkinJacobian) const
//   {
//     // here we assume that the residual policy already:
//     // - reconstucted the previous states
//     // - the projector has already been updated using the galerkinState
//     //	by the residual policy so we don't need to reupdate it.

//     const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
//     const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
//     fomSystemObj.applyDiscreteTimeJacobian(step, time, dt,
// 					   *(phi_.get().data()),
// 					   *fomApplyJac_.data(),
// 					   *yn.data(),
// 					   *ynm1.data());

//     // apply projection
//     projector_.get().apply(fomApplyJac_, galerkinJacobian);
//   }

//   // we have here n = 2 prev rom states
//   template<
//     typename galerkin_state_t,
//     typename galerkin_prev_states_t,
//     typename fom_system_t,
//     typename scalar_t
//     >
//   mpl::enable_if_t< galerkin_prev_states_t::size()==2 >
//   compute_impl(const galerkin_state_t	& galerkinState,
// 	       const galerkin_prev_states_t & galerkinPrevStates,
// 	       const fom_system_t       & fomSystemObj,
// 	       const scalar_t		& time,
// 	       const scalar_t		& dt,
// 	       const time_step_type	& step,
// 	       galerkin_jacobian_type	& galerkinJacobian) const
//   {
//     // here we assume that the residual policy already:
//     // - reconstucted the previous states
//     // - the projector has already been updated using the galerkinState
//     //	by the residual policy so we don't need to reupdate it.

//     const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
//     const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
//     const auto & ynm2 = fomStatesMngr_.get().fomStatePrevStepCRef();
//     fomSystemObj.applyDiscreteTimeJacobian(step, time, dt,
// 					   *(phi_.get().data()),
// 					   *fomApplyJac_.data(),
// 					   *yn.data(),
// 					   *ynm1.data());

//     // apply projection
//     projector_.get().apply(fomApplyJac_, galerkinJacobian);
//   }

// private:
//   const std::size_t romSize_ = {};
//   std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
//   std::reference_wrapper<const decoder_type> decoderObj_ = {};
//   std::reference_wrapper<const typename decoder_type::jacobian_type> phi_;
//   mutable fom_apply_jacobian_ret_type fomApplyJac_;
//   std::reference_wrapper<const projector_t> projector_;
// };

}}}}//end namespace
#endif  // ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_APPLY_JACOBIAN_POLICY_HPP_
