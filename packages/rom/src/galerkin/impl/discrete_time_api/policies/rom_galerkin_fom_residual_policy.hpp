/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_fom_residual_policy.hpp
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

#ifndef ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_RESIDUAL_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <typename fom_states_manager_t, typename fom_residual_type>
class FomResidualPolicyDiscreteTimeApi
{
public:
  using data_type = fom_residual_type;

private:
  // storedStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable ::pressio::ode::types::step_t storedStep_ = {};

  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  mutable fom_residual_type fomResidual_;

public:
  FomResidualPolicyDiscreteTimeApi() = delete;
  FomResidualPolicyDiscreteTimeApi(const FomResidualPolicyDiscreteTimeApi &) = default;
  FomResidualPolicyDiscreteTimeApi & operator=(const FomResidualPolicyDiscreteTimeApi &) = delete;
  FomResidualPolicyDiscreteTimeApi(FomResidualPolicyDiscreteTimeApi &&) = default;
  FomResidualPolicyDiscreteTimeApi & operator=(FomResidualPolicyDiscreteTimeApi &&) = delete;
  ~FomResidualPolicyDiscreteTimeApi() = default;

  template<typename fom_system_t>
  FomResidualPolicyDiscreteTimeApi(const fom_system_t & fomSystemObj,
				fom_states_manager_t & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr),
      fomResidual_(fomSystemObj.createDiscreteTimeResidual())
  {}

public:
  const fom_residual_type & get() const{
    return fomResidual_;
  }

  // we have here n = 1 prev rom states
  template <
    typename galerkin_state_t,
    typename galerkin_prev_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  mpl::enable_if_t< galerkin_prev_states_t::size()==1 >
  compute(const galerkin_state_t & galerkinState,
	  const fom_system_t & fomSystemObj,
	  const scalar_t & time,
	  const galerkin_prev_states_t & galerkinPrevStates,
	  const scalar_t & dt,
	  const ::pressio::ode::types::step_t & step) const
  {
    this->doFomStatesReconstruction(galerkinState, galerkinPrevStates, step);

    const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
    const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
    fomSystemObj.discreteTimeResidual(step, time, dt,
				      *fomResidual_.data(),
				      *yn.data(), *ynm1.data());
  }

  // we have here n = 2 prev rom states
  template <
    typename galerkin_state_t,
    typename galerkin_prev_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  mpl::enable_if_t< galerkin_prev_states_t::size()==2 >
  compute(const galerkin_state_t & galerkinState,
	  const fom_system_t & fomSystemObj,
	  const scalar_t & time,
	  const galerkin_prev_states_t & galerkinPrevStates,
	  const scalar_t & dt,
	  const ::pressio::ode::types::step_t & step) const
  {
    this->doFomStatesReconstruction(galerkinState, galerkinPrevStates, step);

    const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
    const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
    const auto & ynm2 = fomStatesMngr_.get().fomStatePrevPrevStepCRef();
    fomSystemObj.discreteTimeResidual(step, time, dt,
				      *fomResidual_.data(),
				      *yn.data(), *ynm1.data(), *ynm2.data());
  }

private:
  template <typename galerkin_state_t, typename galerkin_prev_states_t>
  void doFomStatesReconstruction(const galerkin_state_t & galerkinState,
				 const galerkin_prev_states_t & galerkinPrevStates,
				 const ::pressio::ode::types::step_t & step) const
  {
    /* the currrent FOM has to be recomputed every time regardless
     * of whether the step changes since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times.
     */
    fomStatesMngr_.get().reconstructCurrentFomState(galerkinState);

    /* previous FOM states should only be recomputed when the time step changes
     * not need to reconstruct all the FOM states, we just need to reconstruct
     * state at the previous step (i.e. t-dt) which is stored in galerkinPrevStates[0]
     */
    if (storedStep_ != step){
      fomStatesMngr_.get() << galerkinPrevStates(ode::nMinusOne());
      storedStep_ = step;
    }
  }
};


// template <
//   typename galerkin_residual_type,
//   typename fom_states_manager_t,
//   typename fom_residual_type,
//   typename projector_t
//   >
// class FomResidualPolicyDiscreteTimeApi
// {
// public:
//   using residual_t = galerkin_residual_type;

// public:
//   FomResidualPolicyDiscreteTimeApi() = delete;
//   FomResidualPolicyDiscreteTimeApi(const FomResidualPolicyDiscreteTimeApi &) = default;
//   FomResidualPolicyDiscreteTimeApi & operator=(const FomResidualPolicyDiscreteTimeApi &) = delete;
//   FomResidualPolicyDiscreteTimeApi(FomResidualPolicyDiscreteTimeApi &&) = default;
//   FomResidualPolicyDiscreteTimeApi & operator=(ResidualPolicyDiscreteTimeApi &&) = delete;
//   ~ResidualPolicyDiscreteTimeApi() = default;

//   template<typename fom_system_t>
//   ResidualPolicyDiscreteTimeApi(std::size_t romSize,
// 				const fom_system_t & fomSystemObj,
// 				fom_states_manager_t & fomStatesMngr,
// 				const projector_t & projector)
//     : romSize_(romSize),
//       fomStatesMngr_(fomStatesMngr),
//       fomResidual_(fomSystemObj.createDiscreteTimeResidual()),
//       projector_(projector)
//   {}

// public:
//   template <typename fom_system_t>
//   galerkin_residual_type create(const fom_system_t & fomSystemObj) const
//   {
//     galerkin_residual_type result(romSize_);
//     ::pressio::ops::set_zero(result);
//     return result;
//   }

//   template <
//     typename ode_tag,
//     typename galerkin_state_t,
//     typename galerkin_prev_states_t,
//     typename fom_system_t,
//     typename scalar_t
//     >
//   void compute(const galerkin_state_t & galerkinState,
// 	       const galerkin_prev_states_t & galerkinPrevStates,
// 	       const fom_system_t & fomSystemObj,
// 	       const scalar_t & time,
// 	       const scalar_t & dt,
// 	       const ::pressio::ode::types::step_t & step,
// 	       galerkin_residual_type & galerkinResidual) const
//   {
//     this->compute_impl(galerkinState, galerkinPrevStates, fomSystemObj,
// 		       time, dt, step, galerkinResidual);
//   }

// private:
//   template <typename galerkin_state_t, typename galerkin_prev_states_t>
//   void doFomStatesReconstruction(const galerkin_state_t & galerkinState,
// 				 const galerkin_prev_states_t & galerkinPrevStates,
// 				 const ::pressio::ode::types::step_t & step) const
//   {
//     /* the currrent FOM has to be recomputed every time regardless
//      * of whether the step changes since we might be inside a non-linear solve
//      * where the time step does not change but this residual method
//      * is called multiple times.
//      */
//     fomStatesMngr_.get().reconstructCurrentFomState(galerkinState);

//     /* previous FOM states should only be recomputed when the time step changes
//      * not need to reconstruct all the FOM states, we just need to reconstruct
//      * state at the previous step (i.e. t-dt) which is stored in galerkinPrevStates[0]
//      */
//     if (storedStep_ != step){
//       fomStatesMngr_.get() << galerkinPrevStates.stateAt(ode::nMinusOne());
//       storedStep_ = step;
//     }
//   }

//   // we have here n = 1 prev rom states
//   template <
//     typename galerkin_state_t,
//     typename galerkin_prev_states_t,
//     typename fom_system_t,
//     typename scalar_t
//     >
//   mpl::enable_if_t< galerkin_prev_states_t::size()==1 >
//   compute_impl(const galerkin_state_t & galerkinState,
// 	       const galerkin_prev_states_t & galerkinPrevStates,
// 	       const fom_system_t & fomSystemObj,
// 	       const scalar_t & time,
// 	       const scalar_t & dt,
// 	       const ::pressio::ode::types::step_t & step,
// 	       galerkin_residual_type & galerkinResidual) const
//   {
//     doFomStatesReconstruction(galerkinState, galerkinPrevStates, step);

//     const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
//     const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
//     fomSystemObj.discreteTimeResidual(step, time, dt,
// 				      *fomResidual_.data(),
// 				      *yn.data(), *ynm1.data());

//     // apply projection
//     projector_.get().apply(fomResidual_, galerkinResidual);
//   }

//   // we have here n = 2 prev rom states
//   template <
//     typename galerkin_state_t,
//     typename galerkin_prev_states_t,
//     typename fom_system_t,
//     typename scalar_t
//     >
//   mpl::enable_if_t< galerkin_prev_states_t::size()==2 >
//   compute_impl(const galerkin_state_t & galerkinState,
// 	       const galerkin_prev_states_t & galerkinPrevStates,
// 	       const fom_system_t & fomSystemObj,
// 	       const scalar_t & time,
// 	       const scalar_t & dt,
// 	       const ::pressio::ode::types::step_t & step,
// 	       galerkin_residual_type & galerkinResidual) const
//   {
//     doFomStatesReconstruction(galerkinState, galerkinPrevStates, step);

//     const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
//     const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
//     const auto & ynm2 = fomStatesMngr_.get().fomStatePrevPrevStepCRef();
//     fomSystemObj.discreteTimeResidual(step, time, dt,
// 				      *fomResidual_.data(),
// 				      *yn.data(), *ynm1.data(), *ynm2.data());

//     // apply projection
//     projector_.get().apply(fomResidual_, galerkinResidual);
//   }

// private:
//   const std::size_t romSize_ = {};

//   // storedStep is used to keep track of which step we are doing.
//   // This is used to decide whether we need to update/recompute the previous
//   // FOM states or not. Since it does not make sense to recompute previous
//   // FOM states if we are not in a new time step.
//   mutable ::pressio::ode::types::step_t storedStep_ = {};

//   std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
//   mutable fom_residual_type fomResidual_;
//   std::reference_wrapper<const projector_t> projector_;
// };

}}}}//end namespace
#endif  // ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_RESIDUAL_POLICY_HPP_
