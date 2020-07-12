/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_jacobian_policy_continuous_time_api.hpp
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

#ifndef ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_VELOCITY_api_HPP_
#define ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_VELOCITY_api_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template<
  typename fom_states_manager_t,
  typename apply_jac_return_type,
  typename decoder_type,
  typename ud_ops_type
  >
class JacobianPolicyContinuousTimeApi
{

public:
  // static constexpr bool isResidualPolicy_ = false;
  using apply_jac_return_t = apply_jac_return_type;
  using ud_ops_t = ud_ops_type;

public:
  JacobianPolicyContinuousTimeApi() = delete;
  ~JacobianPolicyContinuousTimeApi() = default;

  /* for constructing this we need to deal with a few cases
   * 1. void ops
   * 2. non-void ops
   */

  // 1. void ops
  template <
    typename _apply_jac_return_type = apply_jac_return_type,
    typename _ud_ops = ud_ops_type,
    ::pressio::mpl::enable_if_t<std::is_void<_ud_ops>::value, int > =0
    >
  JacobianPolicyContinuousTimeApi(fom_states_manager_t & fomStatesMngr,
			    const _apply_jac_return_type & applyJacObj,
			    const decoder_type & decoder)
    : JJ_(applyJacObj), fomStatesMngr_(fomStatesMngr),  decoderObj_(decoder){}

  // 2. non-void ops
  template <
    typename _apply_jac_return_type = apply_jac_return_type,
    typename _ud_ops = ud_ops_type,
    ::pressio::mpl::enable_if_t<!std::is_void<_ud_ops>::value, int > =0
    >
  JacobianPolicyContinuousTimeApi(fom_states_manager_t & fomStatesMngr,
			    const _apply_jac_return_type & applyJacObj,
			    const decoder_type & decoder,
			    const _ud_ops & udOps)
    : JJ_(applyJacObj),
      fomStatesMngr_(fomStatesMngr),
      decoderObj_(decoder),
      udOps_{&udOps}
  {}

public:

  template <typename system_t>
  apply_jac_return_t create(const system_t & app) const
  {
    return JJ_;
  }

  template <
    typename stepper_tag,
    typename prev_states_mgr,
    typename lspg_state_t, 
    typename lspg_jac_t, 
    typename system_t, 
    typename scalar_t
  >
  void compute(const lspg_state_t & romState,
        const prev_states_mgr & prevStatesMgr,
  		  const system_t	& app,
		  const scalar_t & time,
		  const scalar_t & dt,
		  const ::pressio::ode::types::step_t & step,
		  lspg_jac_t & romJac) const
  {
    this->compute_impl<stepper_tag>(romState, romJac, app, time, dt, step);
  }

private:
  template <
    typename stepper_tag,
    typename matrix_t,
    typename scalar_t,
    typename decoder_jac_type,
    typename _ud_ops = ud_ops_type
  >
  ::pressio::mpl::enable_if_t< std::is_void<_ud_ops>::value >
  time_discrete_dispatcher(matrix_t & romJac,
				scalar_t  dt,
				const decoder_jac_type & phi) const
  {
    ::pressio::rom::lspg::impl::unsteady::time_discrete_jacobian<stepper_tag>(romJac, dt, phi);
  }

  template <
    typename stepper_tag,
    typename matrix_t,
    typename scalar_t,
    typename decoder_jac_type,
    typename _ud_ops = ud_ops_type
  >
  ::pressio::mpl::enable_if_t<!std::is_void<_ud_ops>::value >
  time_discrete_dispatcher(matrix_t & romJac,
				scalar_t dt,
				const decoder_jac_type & phi) const
  {
    ::pressio::rom::lspg::impl::unsteady::time_discrete_jacobian<stepper_tag>(romJac, dt, phi, udOps_);
  }


  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename lspg_jac_t,
    typename system_t,
    typename scalar_t
    >
  void compute_impl(const lspg_state_t & romState,
		    lspg_jac_t	     & romJac,
		    const system_t	     & app,
		    const scalar_t   & t,
		    const scalar_t   & dt,
		    const ::pressio::ode::types::step_t & step) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg apply jac");
#endif

    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    // fomStatesMngr_.reconstructCurrentFomState(romState);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom apply jac");
#endif
    const auto & basis = decoderObj_.getReferenceToJacobian();
    ::pressio::rom::queryFomApplyJacobian(app, fomStatesMngr_.getCRefToCurrentFomState(),
						  basis, romJac, t);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom apply jac");
    timer->start("time discrete jacob");
#endif

    this->time_discrete_dispatcher<stepper_tag>(romJac, dt, basis);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("time discrete jacob");
    timer->stop("lspg apply jac");
#endif
  }


protected:
  mutable apply_jac_return_t JJ_ = {};
  fom_states_manager_t & fomStatesMngr_;
  const decoder_type & decoderObj_ = {};

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  typename std::conditional<
    ::pressio::mpl::is_same<ud_ops_type, pybind11::object>::value, 
    ud_ops_type, const ud_ops_type * >::type udOps_ = {};
#else
    const ud_ops_type * udOps_ = {};
#endif

};

}}}}}
#endif
