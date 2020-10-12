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
  typename wls_state_type, typename decoder_t, typename ode_tag, typename wls_hessian_type, typename policy_t
  >
class SystemHessianGradientApi
{
  // in all cases we need types to be wrappers from containers
  static_assert(::pressio::containers::predicates::is_vector_wrapper<wls_state_type>::value,
		"For WLS, the state_type must be a pressio container vector");
  static_assert(::pressio::containers::predicates::is_dense_matrix_wrapper<wls_hessian_type>::value,
		"WLS: hessian_type must be a pressio container dense matrix");

public:
  // aliases for scalar, state, gradient and hessian NEEDED for solver to detect
  using scalar_type		= typename policy_t::scalar_t;
  using state_type		= wls_state_type;
  using gradient_type		= wls_state_type;
  using hessian_type		= wls_hessian_type;

  using fom_state_t		= typename policy_t::fom_state_t;
  using fom_state_reconstr_t	= ::pressio::rom::FomStateReconstructor<scalar_type, fom_state_t, decoder_t>;

  // information on stencil width, time discrete resiudal, time discrete jacobian, etc.
  using time_stencil_t = ::pressio::rom::wls::timeschemes::timescheme_t<ode_tag, fom_state_t>;
  static constexpr auto timeStencilSize_ = time_stencil_t::state_stencil_size_;

private:
  //size of generalized coordinates
  rom_size_t romSize_ = {};

  //number of discrete time instances in a window
  window_size_t numStepsInWindow_ = {};

  // policy for evaluating the hessian and gradient
  std::reference_wrapper<const policy_t> hessianGradientPolicy_;

  // cache the size of the wls problem: romSize*numStepsInWindow
  window_size_t wlsProblemSize_ = romSize_*numStepsInWindow_;
  window_size_t wlsStencilSize_ = romSize_*timeStencilSize_;

  wls_state_type wlsStateIC_;

  // fom state reconstructor
  const fom_state_reconstr_t fomStateReconstructor_;

  // object knowing the time stencil for doing chuncks of steps
  //time_stencil_t timeSchemeObj_;

  window_size_t activeWindowIndex_ = {};

  // global step number
  window_size_t step_s_		= {};
  scalar_type dt_		= {};
  scalar_type windowStartTime_	= {};

public:
  SystemHessianGradientApi(const rom_size_t romSize,
			   const window_size_t numStepsInWindow,
			   const decoder_t & decoderObj,
			   const policy_t & policy,
			   const fom_state_t & fomStateInitCondition,
			   const fom_state_t & fomNominalState,
			   const wls_state_type & wlsStateIcIn)
    : SystemHessianGradientApi(romSize, numStepsInWindow, decoderObj, policy, fomNominalState)
  {
    const auto spanStartIndex = romSize_*(timeStencilSize_-1);
    auto wlsInitialStateNm1   = containers::span(wlsStateIC_, spanStartIndex, romSize_);
    ::pressio::ops::deep_copy(wlsInitialStateNm1, wlsStateIcIn);
  }

  template <
  typename linear_solver_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::mpl::is_detected<::pressio::solvers::predicates::has_matrix_typedef, linear_solver_type>::value and
    !std::is_void<typename linear_solver_type::matrix_type>::value and
    // has a solveAllowMatOverwrite method
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
  SystemHessianGradientApi(const rom_size_t romSize,
			   const window_size_t numStepsInWindow,
			   const decoder_t & decoderObj,
			   const policy_t & policy,
			   const fom_state_t & fomStateInitCondition,
			   const fom_state_t & fomNominalState,
  			   linear_solver_type & linSolverObj)
    : SystemHessianGradientApi(romSize, numStepsInWindow, decoderObj, policy, fomNominalState)
  {
    // Set initial condition based on L^2 projection onto trial space
    // note that wlsStateIC_[-romSize:end] contains nm1, wlsStateIC[-2*romSize:-romSize] contains nm2 entry, etc.
    wls_state_type wlsStateTmp(romSize);
    const auto & decoderJac = decoderObj.jacobianCRef();

    ::pressio::rom::utils::set_gen_coordinates_L2_projection<scalar_type>(linSolverObj, decoderJac,
  									  fomStateInitCondition, fomNominalState,
									  wlsStateTmp);

    const auto spanStartIndex = romSize_*(timeStencilSize_-1);
    auto wlsInitialStateNm1 = containers::span(wlsStateIC_, spanStartIndex, romSize_);
    ::pressio::ops::deep_copy(wlsInitialStateNm1, wlsStateTmp);
  }

private:
  // delegated constructor to simplify the ones above
  SystemHessianGradientApi(const rom_size_t romSize,
			   const window_size_t numStepsInWindow,
			   const decoder_t & decoderObj,
			   const policy_t & policy,
			   const fom_state_t & fomNominalState)
    : romSize_(romSize),
      numStepsInWindow_{numStepsInWindow},
      hessianGradientPolicy_(policy),
      wlsStateIC_(wlsStencilSize_),
      fomStateReconstructor_(fomNominalState, decoderObj)
  {}

public:
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

  void hessianAndGradient(const state_type	      & wls_state,
			  hessian_type		      & hessian,
			  gradient_type		      & gradient,
			  const pressio::Norm & normType,
			  scalar_type		      & rnorm,
        bool recomputeJacobian) const
  {
    if (normType != ::pressio::Norm::L2)
      throw std::runtime_error("cannot call WLS with a norm != L2");

    rnorm = pressio::utils::constants<scalar_type>::zero();
    hessianGradientPolicy_.get()(
			   wls_state,
			   wlsStateIC_,
			   hessian,
			   gradient,
			   fomStateReconstructor_,
			   dt_,
			   numStepsInWindow_,
			   windowStartTime_,
			   step_s_,
			   rnorm);
  }//end computeHessianAndGradient

  void residualNorm(const state_type	      & wls_state,
		    const pressio::Norm & normType,
		    scalar_type		      & rnorm) const
  {
    if (normType != ::pressio::Norm::L2)
      throw std::runtime_error("cannot call WLS with a norm != L2");

    rnorm = pressio::utils::constants<scalar_type>::zero();
    hessianGradientPolicy_.get().computeResidualNorm(
			   wls_state,
			   wlsStateIC_,
			   fomStateReconstructor_,
			   dt_,
			   numStepsInWindow_,
			   windowStartTime_,
			   step_s_,
			   rnorm);
  }

  // method to advance one window. We may want to put this into some type of window stepper class
  // if we want to have more complex stepping
  template <typename solverType>
  void advanceOneWindow(wls_state_type & wlsState,
			solverType & solver,
			const window_size_t & windowIndex,
			scalar_type dt)
  {
    dt_			= dt;
    activeWindowIndex_  = windowIndex;
    windowStartTime_	= windowIndex*dt_*numStepsInWindow_;
    step_s_		= windowIndex*numStepsInWindow_;

    // set initial guess over window (needed to avoid bad initial guesses that yield NaN)
    for (window_size_t i =0; i < numStepsInWindow_; i++){
      auto wlsViewAssign = ::pressio::containers::span(wlsState,i*romSize_,romSize_);
      auto wlsViewCopy = ::pressio::containers::span(wlsStateIC_,(timeStencilSize_ - 1)*romSize_,romSize_);
      ::pressio::ops::deep_copy(wlsViewAssign, wlsViewCopy);
    }

    // solve system
    solver.solve(*this, wlsState);

    // Loop to update the the wlsStateIC vector.
    // If we add multistep explicit methods, need to add things here.
    const auto start = std::max(0,this->timeStencilSize_  - this->numStepsInWindow_);
    for (window_size_t i = 0; i < start; ++i)
    {
      auto wlsTmpState	      = ::pressio::containers::span(wlsStateIC_, i*romSize_,     romSize_);
      const auto wlsTmpState2 = ::pressio::containers::span(wlsStateIC_, (i+1)*romSize_, romSize_);
      ::pressio::ops::deep_copy(wlsTmpState, wlsTmpState2);
    }

    for (window_size_t i = start ; i < timeStencilSize_; ++i)
    {
      auto wlsTmpState  = ::pressio::containers::span(wlsStateIC_, i*romSize_, romSize_);
      const auto wlsTmpState2 = ::pressio::containers::span(wlsState,
							    (numStepsInWindow_-timeStencilSize_+i)*romSize_,
							    romSize_);
      ::pressio::ops::deep_copy(wlsTmpState, wlsTmpState2);
    }
#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    ::pressio::utils::io::print_stdout("\n");
    auto fmt = ::pressio::utils::io::underline();
    ::pressio::utils::io::print_stdout(fmt, "Completed window ", windowIndex, ::pressio::utils::io::reset(), "\n");
#endif
  }// end advanceOneWindow

};

}}}}//end namespace pressio::rom::src::wls::impl
#endif  // ROM_WLS_IMPL_ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_IMPL_HPP_
