/*
//@HEADER
// ************************************************************************
//
// rom_lspg_jacobian_policy.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef ROM_LSPG_JACOBIAN_POLICY_HPP_
#define ROM_LSPG_JACOBIAN_POLICY_HPP_

#include "../../rom_fwd.hpp"
#include "rom_lspg_time_discrete_jacobian.hpp"
#include "../../../../ode/src/implicit/policies/base/ode_jacobian_policy_base.hpp"
#include "../../rom_data_fom_states.hpp"

namespace pressio{ namespace rom{

template<
  typename fom_states_data,
  typename apply_jac_return_type,
  typename fom_apply_jac_policy,
  typename decoder_type,
  typename ud_ops
  >
class LSPGJacobianPolicy
  : public ode::policy::JacobianPolicyBase<
	LSPGJacobianPolicy<fom_states_data,
			   apply_jac_return_type,
			   fom_apply_jac_policy,
			   decoder_type,
			   ud_ops>>,
    protected fom_states_data,
    protected fom_apply_jac_policy
{

protected:
  using this_t = LSPGJacobianPolicy<fom_states_data,
				    apply_jac_return_type,
				    fom_apply_jac_policy,
				    decoder_type,
				    ud_ops>;

  friend ode::policy::JacobianPolicyBase<this_t>;
  using fom_states_data::yFom_;
  using fom_states_data::yFomOld_;

public:
  static constexpr bool isResidualPolicy_ = false;
  using apply_jac_return_t = apply_jac_return_type;

#ifdef HAVE_PYBIND11
  typename std::conditional<
    mpl::is_same<ud_ops, pybind11::object>::value,
    ud_ops,
    const ud_ops *
    >::type udOps_ = {};
#else
    const ud_ops * udOps_ = {};
#endif


public:
  LSPGJacobianPolicy() = delete;
  ~LSPGJacobianPolicy() = default;

  // this cnstr only enabled when udOps is void
  template <
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      std::is_void<_ud_ops>::value
      > * = nullptr
    >
  LSPGJacobianPolicy(const fom_states_data	 & fomStates,
		     const fom_apply_jac_policy  & applyJacFunctor,
		     const apply_jac_return_type & applyJacObj,
		     const decoder_type		 & decoder)
    : fom_states_data(fomStates),
      fom_apply_jac_policy(applyJacFunctor),
      JJ_(applyJacObj),
      decoderObj_(decoder){
    static_assert( std::is_void<_ud_ops>::value, "");
  }

  // this cnstr only enabled when udOps is non-void
  template <
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      !std::is_void<_ud_ops>::value 
#ifdef HAVE_PYBIND11      
      and mpl::not_same<_ud_ops, pybind11::object>::value
#endif      
      > * = nullptr
    >
  LSPGJacobianPolicy(const fom_states_data	 & fomStates,
		     const fom_apply_jac_policy  & applyJacFunctor,
		     const apply_jac_return_type & applyJacObj,
		     const decoder_type		 & decoder,
		     const _ud_ops & udOps)
    : fom_states_data(fomStates),
      fom_apply_jac_policy(applyJacFunctor),
      JJ_(applyJacObj),
      decoderObj_(decoder),
      udOps_{&udOps}{
    static_assert( !std::is_void<_ud_ops>::value, "");
  }

#ifdef HAVE_PYBIND11
  // this cnstr only enabled when udOps is non-void and python
  template <
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      !std::is_void<_ud_ops>::value and
      mpl::is_same<_ud_ops, pybind11::object>::value
      > * = nullptr
    >
  LSPGJacobianPolicy(const fom_states_data	 & fomStates,
		     const fom_apply_jac_policy  & applyJacFunctor,
		     const apply_jac_return_type & applyJacObj,
		     const decoder_type		 & decoder,
		     const _ud_ops & udOps)
    : fom_states_data(fomStates),
      fom_apply_jac_policy(applyJacFunctor),
      JJ_(applyJacObj),
      decoderObj_(decoder),
      udOps_{udOps}
  {}
#endif
  

public:
  template <ode::ImplicitEnum odeMethod,
	    typename lspg_state_t,
	    typename lspg_jac_t,
	    typename app_t,
	    typename scalar_t>
  void operator()(const lspg_state_t & romY,
		  lspg_jac_t	     & romJJ,
  		  const app_t	     & app,
		  scalar_t	     t,
		  scalar_t	     dt) const
  {
    this->compute_impl<odeMethod>(romY, romJJ, app, t, dt);
  }


  template <ode::ImplicitEnum odeMethod,
	    typename lspg_state_t,
	    typename app_t,
	    typename scalar_t>
  apply_jac_return_t operator()(const lspg_state_t & romY,
				const app_t	   & app,
				scalar_t	   t,
				scalar_t	   dt) const
  {
    this->compute_impl<odeMethod>(romY, JJ_, app, t, dt);
    return JJ_;
  }

private:
  template <
    ode::ImplicitEnum odeMethod,
    typename matrix_t,
    typename scalar_t,
    typename decoder_jac_type,
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
	std::is_void<_ud_ops>::value
      > * = nullptr
  >
  void time_discrete_dispatcher(matrix_t & romJJ,
				scalar_t	dt,
				const decoder_jac_type & phi) const{
    rom::impl::time_discrete_jacobian<odeMethod>(romJJ, dt, phi);
  }

  template <
    ode::ImplicitEnum odeMethod,
    typename matrix_t,
    typename scalar_t,
    typename decoder_jac_type,
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      !std::is_void<_ud_ops>::value
      > * = nullptr
  >
  void time_discrete_dispatcher(matrix_t & romJJ,
				scalar_t dt,
				const decoder_jac_type & phi) const{
    rom::impl::time_discrete_jacobian<odeMethod>(romJJ, dt, phi, udOps_);
  }

  template <ode::ImplicitEnum odeMethod,
	    typename lspg_state_t,
	    typename lspg_jac_t,
	    typename app_t,
	    typename scalar_t>
  void compute_impl(const lspg_state_t & romY,
		    lspg_jac_t	     & romJJ,
		    const app_t	     & app,
		    scalar_t	     t,
		    scalar_t	     dt) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg apply jac");
#endif

    // todo: this is not needed if jacobian is called after resiudal
    // because residual takes care of reconstructing the fom state
    //    timer->start("reconstruct fom state");
    fom_states_data::template reconstructCurrentFomState(romY);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom apply jac");
#endif
    const auto & basis = decoderObj_.getReferenceToJacobian();
    fom_apply_jac_policy::evaluate(app, yFom_, basis, romJJ, t);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("fom apply jac");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("time discrete jacob");
#endif
    this->time_discrete_dispatcher<odeMethod>(romJJ, dt, basis);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time discrete jacob");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg apply jac");
#endif
  }

protected:
  mutable apply_jac_return_t JJ_	= {};
  const decoder_type & decoderObj_	= {};

};

}}//end namespace pressio::rom
#endif
