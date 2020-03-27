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

#ifndef ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_IMPL_HPP_
#define ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_IMPL_HPP_

namespace pressio{ namespace rom{ namespace wls{ namespace impl{

template<
  typename fom_type,
  typename wls_state_type,
  typename decoder_t,
  typename ode_tag,
  typename hessian_t,
  typename hessian_matrix_structure_tag
  >
class SystemHessianGradientApi{

public:
  // aliases for scalar, state, gradient and hessian needed for solver to detect
  // gradient is the same as state for now
  using scalar_type		= typename fom_type::scalar_type;
  using state_type		= wls_state_type;
  using gradient_type		= wls_state_type;
  using hessian_type		= hessian_t;

public:
  using fom_native_state_t	= typename fom_type::state_type;
  using fom_state_t		= ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_reconstr_t	= ::pressio::rom::FomStateReconstructor<scalar_type, fom_state_t, decoder_t>;
  using decoder_jac_t		= typename decoder_t::jacobian_type;

  // policy type (here policy knows how to compute hessian and gradient)
  using hessian_gradient_pol_t	= ::pressio::rom::wls::HessianGradientSequentialPolicy<fom_type,decoder_t,hessian_matrix_structure_tag>;

  // information on stencil width, time discrete resiudal, time discrete jacobian, etc.
  using time_stencil_t = ::pressio::rom::wls::timeschemes::timescheme_t<ode_tag, fom_state_t, wls_state_type>;
  static constexpr auto timeStencilSize_ = time_stencil_t::state_stencil_size_;

  // in all cases we need types to be wrappers from containers
  static_assert(::pressio::containers::meta::is_vector_wrapper<wls_state_type>::value,
		"For WLS, the state_type must be a pressio container vector");
  static_assert(::pressio::containers::meta::is_matrix_wrapper<hessian_type>::value,
		"WLS: hessian_type must be a pressio container matrix");

private:
  const fom_type & appObj_;
  const fom_state_reconstr_t  fomStateReconstructorObj_;

  //size of generalized coordinates
  int romSize_			= {};
  // object knowing the time stencil for doing chuncks of step
  time_stencil_t timeSchemeObj_;

  //number of discrete time instances in a window
  int numStepsInWindow_		= {};

  // policy for evaluating the hessian and gradient
  const hessian_gradient_pol_t hessian_gradient_polObj_;

  int activeWindowIndex_	= {};
  // global step number
  int step_s_			= {};
  scalar_type dt_		= {};
  scalar_type windowStartTime_	= {};

  // cache the size of the wls problem: romSize*numStepsInWindow
  int wlsProblemSize_		= romSize_*numStepsInWindow_;
  int wlsStencilSize_           = romSize_*timeStencilSize_;
  // I keep this here since you pass it to the policy. Originally this was owened
  // by the timeSchemeObj but it should be owened here
  wls_state_type wlsStateIC_;

public:

  SystemHessianGradientApi(const fom_type & appObj,
			   const fom_state_t & yFOM_IC,
			   const fom_state_t & yFOM_Ref,
			   const decoder_t & decoderObj,
			   const int numStepsInWindow,
			   const int romSize,
			   const wls_state_type & wlsStateIcIn)
    : appObj_(appObj),
      fomStateReconstructorObj_(yFOM_Ref, decoderObj),
      romSize_(romSize),
      timeSchemeObj_(romSize_, yFOM_IC),
      numStepsInWindow_{numStepsInWindow},
      hessian_gradient_polObj_( appObj, yFOM_IC, numStepsInWindow, timeStencilSize_, decoderObj,romSize),
      wlsStateIC_(wlsStencilSize_)
  {
    const auto spanStartIndex = romSize_*(timeStencilSize_-1);
    auto wlsInitialStateNm1 = containers::span(wlsStateIC_, spanStartIndex, romSize_);
    ::pressio::ops::deep_copy(wlsInitialStateNm1, wlsStateIcIn);
  }


  template <
    typename linear_solver_type,
    ::pressio::mpl::enable_if_t<
      ::pressio::mpl::is_detected<::pressio::solvers::meta::has_matrix_typedef, linear_solver_type>::value and
      !std::is_void<typename linear_solver_type::matrix_type>::value and
      ::pressio::mpl::publicly_inherits_from<
	linear_solver_type,
	::pressio::solvers::LinearBase<
	  typename linear_solver_type::matrix_type, linear_solver_type
	  >
	>::value
      > * = nullptr
  >
  SystemHessianGradientApi(const fom_type & appObj,
			   const fom_state_t & yFOM_IC,
			   const fom_state_t & yFOM_Ref,
			   const decoder_t & decoderObj,
			   const int numStepsInWindow,
			   const int romSize,
			   linear_solver_type & linSolverObj)
    : appObj_(appObj),
      fomStateReconstructorObj_(yFOM_Ref, decoderObj),
      romSize_(romSize),
      timeSchemeObj_(romSize_, yFOM_IC),
      numStepsInWindow_{numStepsInWindow},
      hessian_gradient_polObj_( appObj, yFOM_IC, numStepsInWindow, timeStencilSize_, decoderObj,romSize),
      wlsStateIC_(wlsStencilSize_)
  {
    // Set initial condition based on L^2 projection onto trial space
    // note that wlsStateIC_[-romSize:end] contains nm1, wlsStateIC[-2*romSize:-romSize] contains nm2 entry, etc.
    wls_state_type wlsStateTmp(romSize);
    const auto & decoderJac = decoderObj.getReferenceToJacobian();

    ::pressio::rom::utils::set_gen_coordinates_L2_projection<scalar_type>(linSolverObj, decoderJac,
									  yFOM_IC, yFOM_Ref, wlsStateTmp);

    const auto spanStartIndex = romSize_*(timeStencilSize_-1);
    auto wlsInitialStateNm1 = containers::span(wlsStateIC_, spanStartIndex, romSize_);
    ::pressio::ops::deep_copy(wlsInitialStateNm1, wlsStateTmp);
  }

  hessian_type createHessianObject(const wls_state_type & stateIn) const{
    hessian_type H(wlsProblemSize_, wlsProblemSize_);
    return H;
  }

  gradient_type createGradientObject(const wls_state_type & stateIn) const{
    gradient_type g(wlsProblemSize_);
    return g;
  }

  void computeHessianAndGradient(const state_type	      & wls_state,
                                 hessian_type		      & hessian,
                                 gradient_type		      & gradient,
                                 const pressio::solvers::Norm & normType,
                                 scalar_type		      & rnorm) const
  {
    rnorm = pressio::utils::constants::zero<scalar_type>();
    hessian_gradient_polObj_(timeSchemeObj_,
                             wls_state,
                             wlsStateIC_,
                             hessian,
                             gradient,
                             fomStateReconstructorObj_,
                             dt_,
                             numStepsInWindow_,
                             windowStartTime_,
                             step_s_,
                             rnorm);
  }//end computeHessianAndGradient

  // method to advance one window. We may want to put this into some type of window stepper class
  // if we want to have more complex stepping
  template <typename solverType>
  void advanceOneWindow(wls_state_type & wlsState,
			solverType & solver,
			const int & windowIndex,
			scalar_type dt)
  {
    #ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      const int procRank = Teuchos::GlobalMPISession::getRank();
      Teuchos::oblackholestream blackhole;
      std::ostream &out = (procRank == 0 ? std::cout : blackhole);
      Teuchos::RCP<Teuchos::StackedTimer> timer_ = Teuchos::TimeMonitor::getStackedTimer();
    #endif

    dt_			= dt; //set time step
    activeWindowIndex_  = windowIndex; //set window number
    windowStartTime_	= windowIndex*dt_*numStepsInWindow_;  //set starting time
    step_s_		= windowIndex*numStepsInWindow_;  //set step number

    solver.solve(*this, wlsState); //solve system

    // Loop to update the the wlsStateIC vector.
    // If we add multistep explicit methods, need to add things here.
    const int start = std::max(0,this->timeStencilSize_  - this->numStepsInWindow_);
    for (int i = 0; i < start; i++)
    {
      auto wlsTmpState	      = ::pressio::containers::span(wlsStateIC_, i*romSize_,     romSize_);
      const auto wlsTmpState2 = ::pressio::containers::span(wlsStateIC_, (i+1)*romSize_, romSize_);
      ::pressio::ops::deep_copy(wlsTmpState, wlsTmpState2);
    }

    for (int i = start ; i < timeStencilSize_; i++)
    {
      auto wlsTmpState  = ::pressio::containers::span(wlsStateIC_, i*romSize_, romSize_);
      const auto wlsTmpState2 = ::pressio::containers::span(wlsState, (numStepsInWindow_ - timeStencilSize_+i)*romSize_, romSize_);
      ::pressio::ops::deep_copy(wlsTmpState, wlsTmpState2);
    }
    std::cout << " Window " << windowIndex << " completed " << std::endl;

    #ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer_->report(out);
    #endif

  }// end advanceOneWindow

};

}}}}//end namespace pressio::rom::src::wls::impl
#endif
