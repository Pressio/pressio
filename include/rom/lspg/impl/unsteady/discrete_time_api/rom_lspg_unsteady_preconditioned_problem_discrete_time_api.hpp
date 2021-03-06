/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_preconditioned_problem_discrete_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_ROM_LSPG_UNSTEADY_PRECONDITIONED_PROBLEM_DISCRETE_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_ROM_LSPG_UNSTEADY_PRECONDITIONED_PROBLEM_DISCRETE_TIME_API_HPP_


namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template<typename ...Args>
class PreconditionedProblemDiscreteTimeApi
{
public:
  using this_t = PreconditionedProblemDiscreteTimeApi<Args...>;
  using traits = ::pressio::rom::details::traits<this_t>;

  using fom_system_t		= typename traits::fom_system_t;
  using scalar_t		= typename traits::scalar_t;
  using fom_native_state_t	= typename traits::fom_native_state_t;
  using fom_native_residual_t	= typename traits::fom_native_residual_t;
  using fom_state_t		= typename traits::fom_state_t;
  using decoder_t		= typename traits::decoder_t;
  using fom_state_reconstr_t	= typename traits::fom_state_reconstr_t;
  using fom_states_manager_t	= typename traits::fom_states_manager_t;
  using preconditioner_t  = typename traits::preconditioner_t;
  using ud_ops_t		= typename traits::ud_ops_t;
  using lspg_state_t		= typename traits::lspg_state_t;
  using residual_policy_t	= typename traits::residual_policy_t;
  using jacobian_policy_t	= typename traits::jacobian_policy_t;
  using stepper_t		= typename traits::stepper_t;

private:
  using At = ::pressio::rom::impl::FomObjMixin<fom_system_t>;
  using Bt = ::pressio::rom::impl::FomStatesMngrMixin<At, ud_ops_t, fom_state_t,
				fom_state_reconstr_t, fom_states_manager_t>;
  using Ct = PrecondPoliciesMixin<Bt, ud_ops_t, residual_policy_t, jacobian_policy_t>;
  using mem_t = ::pressio::rom::impl::ImplicitStepperMixin<Ct, void, stepper_t>;
  mem_t members_;

public:
  stepper_t & stepperRef(){ return members_.stepperObj_; }

  const fom_native_state_t & currentFomStateCRef() const{
    return *(members_.fomStatesMngr_(::pressio::ode::nPlusOne()).data());
  }

  const fom_state_reconstr_t & fomStateReconstructorCRef() const{
    return members_.fomStateReconstructor_;
  }

public:
  PreconditionedProblemDiscreteTimeApi() = delete;
  PreconditionedProblemDiscreteTimeApi(const PreconditionedProblemDiscreteTimeApi &) = default;
  PreconditionedProblemDiscreteTimeApi & operator=(const PreconditionedProblemDiscreteTimeApi &) = delete;
  PreconditionedProblemDiscreteTimeApi(PreconditionedProblemDiscreteTimeApi &&) = default;
  PreconditionedProblemDiscreteTimeApi & operator=(PreconditionedProblemDiscreteTimeApi &&) = delete;
  ~PreconditionedProblemDiscreteTimeApi() = default;

  /* ud_ops_t == void */
  template<
    typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t< std::is_void<_ud_ops_t>::value, int > = 0
    >
  PreconditionedProblemDiscreteTimeApi(const fom_system_t & fomSystemObj,
				       decoder_t & decoder,
				       const lspg_state_t & romStateIn,
				       const fom_native_state_t & fomNominalStateNative,
				       const preconditioner_t & preconditionerObj)
    : members_(romStateIn, fomSystemObj, decoder,
	       fomNominalStateNative, preconditionerObj)
    {}
};

}}}}}//end namespace pressio::rom::lspg::unsteady::impl
#endif  // ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_ROM_LSPG_UNSTEADY_PRECONDITIONED_PROBLEM_DISCRETE_TIME_API_HPP_
