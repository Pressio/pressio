/*
//@HEADER
// ************************************************************************
//
// rom_wls_hessian_gradient_system_api_impl.hpp
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

#ifndef ROM_WLS_IMPL_ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_IMPL_HPP_
#define ROM_WLS_IMPL_ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_IMPL_HPP_

namespace pressio{ namespace rom{ namespace wls{ namespace impl{

template<
  typename wls_state_type,
  typename decoder_t,
  typename wls_hessian_type,
  typename policy_t
>
class SystemHessianGradientApi
{
  static_assert
  (::pressio::containers::predicates::is_wrapper<wls_state_type>::value and
   wls_state_type::traits::rank == 1,
   "For WLS, the state_type must be a rank-1 pressio container");

  static_assert
  (::pressio::containers::predicates::is_wrapper<wls_hessian_type>::value and
   wls_hessian_type::traits::rank == 2,
   "WLS: hessian_type must be a rank-2 pressio container");

public:
  using scalar_type		= typename policy_t::scalar_t;
  using state_type		= wls_state_type;
  using gradient_type		= wls_state_type;
  using hessian_type		= wls_hessian_type;
  using fom_state_t		= typename decoder_t::fom_state_type;
  using fom_native_state_t      = typename fom_state_t::traits::wrapped_t;
  using fom_state_reconstr_t	= ::pressio::rom::FomStateReconstructor<scalar_type, fom_state_t, decoder_t>;
  using wls_native_state_t	= typename wls_state_type::traits::wrapped_t;

public:
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  // constructor needed for pressio4py
  SystemHessianGradientApi(const decoder_t & decoderObj,
			   const policy_t & policy,
			   const fom_native_state_t & fomStateInitCondition,
			   const fom_native_state_t & fomNominalState,
			   const wls_native_state_t & wlsStateIcIn)
    : wlsProblemSize_(policy.romSize()*policy.numStepsInWindow()),
      hessianGradientPolicy_(policy),
      wlsStateIC_(policy_t::timeStencilSize_*policy.romSize()),
      fomStateReconstructor_(fomNominalState, decoderObj)
  {
    const auto spanStartIndex = policy.romSize()*(policy_t::timeStencilSize_-1);
    auto wlsInitialStateNm1 = containers::span(wlsStateIC_, spanStartIndex, policy.romSize());
    ::pressio::ops::deep_copy(wlsInitialStateNm1,
			      wls_state_type(wlsStateIcIn, ::pressio::view()));
  }

#else

  SystemHessianGradientApi(const decoder_t & decoderObj,
			   const policy_t & policy,
			   const fom_state_t & fomStateInitCondition,
			   const fom_state_t & fomNominalState,
			   const wls_state_type & wlsStateIcIn)
    : SystemHessianGradientApi(decoderObj, policy, fomNominalState)
  {
    const auto spanStartIndex = policy.romSize()*(policy_t::timeStencilSize_-1);
    auto wlsInitialStateNm1   = containers::span(wlsStateIC_, spanStartIndex, policy.romSize());
    ::pressio::ops::deep_copy(wlsInitialStateNm1, wlsStateIcIn);
  }

  template <
  typename linear_solver_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::mpl::is_detected<
      ::pressio::solvers::predicates::has_matrix_typedef, linear_solver_type
      >::value and
    !std::is_void<typename linear_solver_type::matrix_type>::value and
    std::is_void<
      decltype
      (
       std::declval<linear_solver_type>().solveAllowMatOverwrite
       (
        std::declval<typename linear_solver_type::matrix_type &>(),
        std::declval<wls_state_type const &>(),
        std::declval<wls_state_type &>()
        )
       )
      >::value,
      int > = 0
  >
  SystemHessianGradientApi(const decoder_t & decoderObj,
			   const policy_t & policy,
			   const fom_state_t & fomStateInitCondition,
			   const fom_state_t & fomNominalState,
  			   linear_solver_type & linSolverObj)
    : SystemHessianGradientApi(decoderObj, policy, fomNominalState)
  {
    // Set initial condition based on L^2 projection onto trial space
    // note that wlsStateIC_[-romSize:end] contains nm1, wlsStateIC[-2*romSize:-romSize] contains nm2 entry, etc.
    wls_state_type wlsStateTmp(policy.romSize());
    const auto & decoderJac = decoderObj.jacobianCRef();

    ::pressio::rom::utils::set_gen_coordinates_L2_projection<scalar_type>(linSolverObj,
									  decoderJac,
									  fomStateInitCondition,
									  fomNominalState,
									  wlsStateTmp);

    const auto spanStartIndex = policy.romSize()*(policy_t::timeStencilSize_-1);
    auto wlsInitialStateNm1 = containers::span(wlsStateIC_, spanStartIndex, policy.romSize());
    ::pressio::ops::deep_copy(wlsInitialStateNm1, wlsStateTmp);
  }
#endif

private:
  // delegated constructor to simplify the ones above
  SystemHessianGradientApi(const decoder_t & decoderObj,
			   const policy_t & policy,
			   const fom_state_t & fomNominalState)
    : wlsProblemSize_(policy.romSize()*policy.numStepsInWindow()),
      hessianGradientPolicy_(policy),
      wlsStateIC_(policy_t::timeStencilSize_*policy.romSize()),
      fomStateReconstructor_(fomNominalState, decoderObj)
  {}

public:
  const ::pressio::rom::wls::window_size_t timeStencilSize() const{
    return policy_t::timeStencilSize_;
  }

  const fom_state_reconstr_t & fomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

  hessian_type createHessian() const{
    hessian_type H(wlsProblemSize_, wlsProblemSize_);
    return H;
  }

  gradient_type createGradient() const{
    gradient_type g(wlsProblemSize_);
    return g;
  }

  void hessianAndGradient(const state_type & wls_state,
			  hessian_type & hessian,
			  gradient_type & gradient,
			  const pressio::Norm & normType,
			  scalar_type & rnorm,
			  bool recomputeJacobian) const
  {
    if (normType != ::pressio::Norm::L2){
      throw std::runtime_error("cannot call WLS with a norm != L2");
    }

    rnorm = pressio::utils::constants<scalar_type>::zero();
    hessianGradientPolicy_(wls_state,
			   wlsStateIC_,
			   hessian,
			   gradient,
			   fomStateReconstructor_,
			   dt_,
			   windowStartTime_,
			   step_s_,
			   rnorm);
  }

  void residualNorm(const state_type & wls_state,
		    const pressio::Norm & normType,
		    scalar_type	& rnorm) const
  {
    if (normType != ::pressio::Norm::L2){
      throw std::runtime_error("cannot call WLS with a norm != L2");
    }

    rnorm = pressio::utils::constants<scalar_type>::zero();
    hessianGradientPolicy_.get().computeResidualNorm
      (
       wls_state,
       wlsStateIC_,
       fomStateReconstructor_,
       dt_,
       windowStartTime_,
       step_s_,
       rnorm);
  }

  void setTimeStepSize(scalar_type dtIn){
    dt_ = dtIn;
  }

  void setWindowStartTime(scalar_type value){
    windowStartTime_ = value;
  }

  void setStepS(window_size_t value){
    step_s_ = value;
  }

  window_size_t numStepsInWindow() const{
    return hessianGradientPolicy_.get().numStepsInWindow();
  }

  rom_size_t romSize() const{
    return hessianGradientPolicy_.get().romSize();
  }

  wls_state_type & wlsStateInitConditionRef(){
    return wlsStateIC_;
  }

private:
  // cache the size of the wls problem: romSize*numStepsInWindow
  window_size_t wlsProblemSize_ = {};

  // policy for evaluating the hessian and gradient
  std::reference_wrapper<const policy_t> hessianGradientPolicy_;

  wls_state_type wlsStateIC_;

  // fom state reconstructor
  const fom_state_reconstr_t fomStateReconstructor_;

  // global step number
  window_size_t step_s_		= {};
  scalar_type dt_		= {};
  scalar_type windowStartTime_	= {};
};

}}}}//end namespace pressio::rom::src::wls::impl
#endif  // ROM_WLS_IMPL_ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_IMPL_HPP_
