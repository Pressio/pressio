/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_fom_apply_jacobian_policy.hpp
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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_FOM_APPLY_JACOBIAN_POLICY_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_FOM_APPLY_JACOBIAN_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <class fom_states_manager_t, class fom_apply_jac_type, class decoder_type>
class FomApplyJacobianPolicy
{
public:
  using data_type = fom_apply_jac_type;

private:
  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  std::reference_wrapper<const typename decoder_type::jacobian_type> phi_;
  mutable fom_apply_jac_type fomApplyJac_ = {};

public:
  FomApplyJacobianPolicy() = delete;
  FomApplyJacobianPolicy(const FomApplyJacobianPolicy &) = default;
  FomApplyJacobianPolicy & operator=(const FomApplyJacobianPolicy &) = delete;
  FomApplyJacobianPolicy(FomApplyJacobianPolicy &&) = default;
  FomApplyJacobianPolicy & operator=(FomApplyJacobianPolicy &&) = delete;
  ~FomApplyJacobianPolicy() = default;

  template<typename fom_system_t>
  FomApplyJacobianPolicy(const fom_system_t & fomSystemObj,
			 fom_states_manager_t & fomStatesMngr,
			 const decoder_type & decoder)
    : fomStatesMngr_(fomStatesMngr),
      phi_(decoder.jacobianCRef()),
      fomApplyJac_(fomSystemObj.createApplyJacobianResult(*phi_.get().data()))
  {}

public:
  template<class galerkin_state_t, typename fom_system_t, typename scalar_t>
  void compute(const galerkin_state_t & galerkinState,
	       const fom_system_t  & fomSystemObj,
	       const scalar_t & evaluationTime) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("galerkin apply fom jacobian");
#endif
    /* here we assume that current FOM state has been reconstructd
       by the residual policy. So we do not recompute the FOM state.
       Maybe we should find a way to ensure this is the case.
     */
    const auto & fomState = fomStatesMngr_.get().fomStateAt(::pressio::ode::nPlusOne());

    // call applyJacobian on the fom object
    fomSystemObj.applyJacobian(*fomState.data(),
			       *phi_.get().data(),
			       evaluationTime,
			       *fomApplyJac_.data());

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("galerkin apply fom jacobian");
#endif
  }

  const fom_apply_jac_type & get() const{ return fomApplyJac_; }
};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_FOM_APPLY_JACOBIAN_POLICY_HPP_
