/*
//@HEADER
// ************************************************************************
//
// rom_wls_hessian_and_gradient_sequential_policy_impl.hpp
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

#ifndef ROM_WLS_IMPL_POLICIES_ROM_WLS_HESSIAN_AND_GRADIENT_SEQUENTIAL_POLICY_IMPL_HPP_
#define ROM_WLS_IMPL_POLICIES_ROM_WLS_HESSIAN_AND_GRADIENT_SEQUENTIAL_POLICY_IMPL_HPP_

/*
Class used for computing the hessian and gradients in the WLS system.

The hessian_gradient policy is responsible for assembling Jw^T Jw and J^T rw,
where Jw is the reduced windowed Jacobian and rw the residual.
The policy initializes in memory:
- wlsJacs: a vector container of n_s -1 local Jacobians, J, where n_s is
  the width of the time stencil (e.g., n_s = 3 for BDF2)
- residual: this is the residual vector for the FOM
- yFOM_current: this is a working variable for the fom state.
*/

namespace pressio{ namespace rom{ namespace wls{ namespace impl{

template<
  typename fom_system_type,
  typename decoder_t,
  typename ode_tag,
  typename hessian_matrix_structure_tag,
  typename preconditioner_t,
  typename jacobians_container_t
  >
class HessianGradientSequentialPolicy
{

public:
  using scalar_t           = typename fom_system_type::scalar_type;
  using fom_state_t	   = typename decoder_t::fom_state_type;
  using fom_native_state_t = typename fom_system_type::state_type;
  using decoder_jac_t      = typename decoder_t::jacobian_type;
  using residual_t         = fom_state_t;
  using time_stencil_t	   = ::pressio::rom::wls::timeschemes::timescheme_t<ode_tag, fom_state_t>;
  static constexpr auto timeStencilSize_ = time_stencil_t::state_stencil_size_;

public:

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  // Constructor for continuous-time API for pressio4py
  template<
  typename U = fom_system_type,
  mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<U>::value
    > * = nullptr
  >
  HessianGradientSequentialPolicy(rom_size_t romSize,
				  window_size_t numStepsInWindow,
				  const decoder_t & decoderObj,
				  pybind11::object fomSystemObj,
				  const fom_native_state_t & fomState,
				  const window_size_t jacobianUpdateFrequency = 1)
    : romSize_(romSize),
      numStepsInWindow_(numStepsInWindow),
      jacStencilSize_(std::min(timeStencilSize_+1, numStepsInWindow)),
      jacobianUpdateFrequency_(jacobianUpdateFrequency),
      fomSystemObj_(fomSystemObj),
      phi_(decoderObj.jacobianCRef()),
      fomStateCurrent_(fomState),
      residual_(fomSystemObj_.fomCRef().createVelocity()),
      J_( fomSystemObj_.fomCRef().createApplyJacobianResult(*(phi_.get().data())) ),
      // construct wls Jacobians from jacobian of the decoder: we might need to change this later
      jacobians_(timeStencilSize_, numStepsInWindow, J_),
      timeSchemeObj_(romSize_, fomState)
  {
    checkJacobianUpdatingFrequency(jacobianUpdateFrequency);
  }

#else

  // Constructor for continuous-time API
  template<
    typename U = fom_system_type,
    mpl::enable_if_t<
      ::pressio::rom::constraints::most_likely_continuous_time_system<U>::value
      > * = nullptr
    >
  HessianGradientSequentialPolicy(rom_size_t romSize,
				  window_size_t numStepsInWindow,
				  const decoder_t & decoderObj,
				  const fom_system_type & fomSystemObj,
				  const fom_state_t & fomState,
				  const window_size_t jacobianUpdateFrequency = 1)
    : romSize_(romSize),
      numStepsInWindow_(numStepsInWindow),
      jacStencilSize_(std::min(timeStencilSize_+1, numStepsInWindow)),
      jacobianUpdateFrequency_(jacobianUpdateFrequency),
      fomSystemObj_(fomSystemObj),
      phi_(decoderObj.jacobianCRef()),
      fomStateCurrent_(fomState),
      residual_(fomSystemObj.createVelocity()),
      J_( fomSystemObj_.fomCRef().createApplyJacobianResult(*(phi_.get().data())) ),
      // construct wls Jacobians from jacobian of the decoder: we might need to change this later
      jacobians_( timeStencilSize_, numStepsInWindow ,J_),
      timeSchemeObj_(romSize_, fomState)
  {
    checkJacobianUpdatingFrequency(jacobianUpdateFrequency);
  }

  // constructor for discrete-time API
  template<
    typename U = fom_system_type,
    mpl::enable_if_t<
      ::pressio::rom::constraints::most_likely_discrete_time_system<U>::value
      > * = nullptr
    >
  HessianGradientSequentialPolicy(rom_size_t romSize,
				  window_size_t numStepsInWindow,
				  const decoder_t & decoderObj,
				  const fom_system_type & fomSystemObj,
				  const fom_state_t & fomState,
				  const window_size_t jacobianUpdateFrequency = 1)
    : romSize_(romSize),
      numStepsInWindow_(numStepsInWindow),
      jacStencilSize_(std::min(timeStencilSize_+1, numStepsInWindow)),
      jacobianUpdateFrequency_(jacobianUpdateFrequency),
      fomSystemObj_(fomSystemObj),
      phi_(decoderObj.jacobianCRef()),
      fomStateCurrent_(fomState),
      residual_( fomSystemObj.createDiscreteTimeResidual() ),
      J_(fomSystemObj.createApplyDiscreteTimeJacobianResult(*(phi_.get().data())) ),
      // construct wls Jacobians from jacobian of the decoder: we might need to change this later
      jacobians_( timeStencilSize_, numStepsInWindow , J_),
      timeSchemeObj_(romSize_, fomState)
  {
    checkJacobianUpdatingFrequency(jacobianUpdateFrequency);
  }
#endif

public:
  rom_size_t romSize() const{
    return romSize_;
  }

  window_size_t numStepsInWindow() const{
    return numStepsInWindow_;
  }

  template <
  typename wls_state_type,
  typename fom_state_reconstr_t,
  typename hess_type,
  typename gradient_type
  >
  void operator()(const wls_state_type  & wlsState,
                  const wls_state_type & wlsStateIC,
                  hess_type & hess,
                  gradient_type & gradient,
                  const fom_state_reconstr_t & fomStateReconstrObj,
                  const scalar_t dt,
                  const scalar_t ts ,
                  const window_size_t step_s,
                  scalar_t & rnorm) const
  {
    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();

    window_size_t stepNumLocal = 0;
    scalar_t t = ts + stepNumLocal*dt;
    window_size_t stepNumGlobal = step_s + stepNumLocal;
    ::pressio::ops::set_zero(hess);
    ::pressio::ops::set_zero(gradient);

    //get access to the state at the first window
    setCurrentFomState(wlsState, 0, fomStateReconstrObj);

    //reconstruct the FOM states from the previous window/ICs
    timeSchemeObj_.updateStatesFirstStep(wlsStateIC, fomStateReconstrObj);

    // ----------------------------------------------------
    //compute the time discrete residual
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("residual");
#endif
    timeSchemeObj_.time_discrete_residual(fomSystemObj_.fomCRef(),
					  fomStateCurrent_, residual_,
					  ts, dt, stepNumGlobal);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("residual");
#endif

    // ----------------------------------------------------
    // preconditioner
    Preconditioner(fomSystemObj_.fomCRef(), fomStateCurrent_, residual_, t);

    //increment the norm
    rNormHelper_ = ::pressio::ops::norm2(residual_);
    rnorm += rNormHelper_*rNormHelper_;

    // ----------------------------------------------------
    // compute jacobian over stencil
    if (innerLoopCounter_ % jacobianUpdateFrequency_ == 0)
    {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("jacobian");
#endif
      computeJacobiansOverStencil(timeSchemeObj_, stepNumGlobal, stepNumLocal , t,dt);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("jacobian");
#endif
    }

    // ----------------------------------------------------
    // compute local block of hessian
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("hessian computation");
#endif

    auto hess_block = ::pressio::containers::subspan
      (hess,
       std::make_pair( stepNumLocal*romSize_,(stepNumLocal+1)*romSize_ ) ,
       std::make_pair( stepNumLocal*romSize_,(stepNumLocal+1)*romSize_) );

    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(), one,
			    jacobians_.localJacobian(stepNumLocal,0),
			    jacobians_.localJacobian(stepNumLocal,0),
			    zero, hess_block);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("hessian computation");
#endif

    // ----------------------------------------------------
    // compute gradient[n*romSize_:(n+1)*romSize] += J^T r
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("gradient computation");
#endif

    auto gradientView = ::pressio::containers::span(gradient, stepNumLocal*romSize_, romSize_);
    ::pressio::ops::product(::pressio::transpose(), one,
			    jacobians_.localJacobian(stepNumLocal,0),
			    residual_,
			    one, gradientView);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("gradient computation");
#endif

    // ----------------------------------------------------
    // loop over local steps in window
    for (window_size_t stepNumLocal = 1; stepNumLocal < numStepsInWindow_; stepNumLocal++)
    {
      updateResidualAndJacobian(timeSchemeObj_, wlsState,
				fomStateReconstrObj,
				stepNumLocal,
				step_s, ts,
				rnorm, gradient, dt);

      const window_size_t sbar = std::min(stepNumLocal+1, jacStencilSize_);
      for (window_size_t i=0; i < sbar; i++)
	{
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
	  timer->start("gradient computation");
#endif
	  auto gradientView = ::pressio::containers::span(gradient,
							  (stepNumLocal-i)*romSize_,
							  romSize_);

	  ::pressio::ops::product(::pressio::transpose(), one,
				  jacobians_.localJacobian( stepNumLocal, i ),
				  residual_, one,
				  gradientView);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
	  timer->stop("gradient computation");
#endif
	}

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("hessian computation");
#endif
      addToHessian(hess, stepNumLocal, sbar, hessian_matrix_structure_tag());
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("hessian computation");
#endif
    }//end loop over stepsInWindow

    innerLoopCounter_ += 1;
    rnorm = std::sqrt(rnorm);

  }//end operator()


  template <
    typename wls_state_type,
    typename fom_state_reconstr_t
    >
  void computeResidualNorm(const wls_state_type  & wlsState,
			   const wls_state_type & wlsStateIC,
			   const fom_state_reconstr_t & fomStateReconstrObj,
			   const scalar_t dt,
			   const scalar_t ts ,
			   const window_size_t step_s,
			   scalar_t & rnorm) const
  {

    window_size_t stepNumLocal = 0;
    scalar_t t = ts + stepNumLocal*dt;
    window_size_t stepNumGlobal = step_s + stepNumLocal;
    //get access to the state at the first window
    setCurrentFomState(wlsState, 0, fomStateReconstrObj);

    //reconstruct the FOM states from the previous window/ICs
    timeSchemeObj_.updateStatesFirstStep(wlsStateIC, fomStateReconstrObj);

    //compute the time discrete residual
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("residual");
#endif

    timeSchemeObj_.time_discrete_residual(fomSystemObj_.fomCRef(),
					  fomStateCurrent_, residual_,
					  ts, dt, stepNumGlobal);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("residual");
#endif
    Preconditioner(fomSystemObj_.fomCRef(), fomStateCurrent_, residual_, t);

    //increment the norm
    rNormHelper_ = ::pressio::ops::norm2(residual_);
    rnorm += rNormHelper_*rNormHelper_;

    for (window_size_t stepNumLocal = 1; stepNumLocal < numStepsInWindow_; stepNumLocal++)
    {
      // reconstruct FOM states
      timeSchemeObj_.updateStatesNStep(fomStateCurrent_);
      setCurrentFomState(wlsState,stepNumLocal,fomStateReconstrObj);

      // Evaluate residual
      t = ts + stepNumLocal*dt;
      window_size_t step;
      step = step_s + stepNumLocal;

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      auto timer = Teuchos::TimeMonitor::getStackedTimer();
      timer->start("residual");
#endif
      timeSchemeObj_.time_discrete_residual(fomSystemObj_.fomCRef(),
					    fomStateCurrent_, residual_,
					    t, dt, step);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("residual");
#endif
      rNormHelper_ = ::pressio::ops::norm2(residual_);
      rnorm += rNormHelper_*rNormHelper_;

    }//end loop over stepsInWindow

    rnorm = std::sqrt(rnorm);

  }//end computeResidualNorm()

private:
  void checkJacobianUpdatingFrequency(const window_size_t jacobianUpdateFrequency)
  {
    using non_frozen_t = ::pressio::rom::wls::NonFrozenJacobiansContainer<typename decoder_t::jacobian_type>;
    constexpr auto nonFrozen = std::is_same<jacobians_container_t, non_frozen_t>::value;

    if (nonFrozen and jacobianUpdateFrequency != 1 )
    {
      PRESSIOLOG_WARN("Warning: Using NonFrozenJacobianContainer with jacobianUpdateFrequency > 1.\
To run with jacobianUpdateFrequency > 1, use FrozenJacobiansContainer\n");
      PRESSIOLOG_WARN("Setting jacobianUpdateFrequency = 1. \n");
      jacobianUpdateFrequency_ = 1;
    }
  }

  template<typename hessian_type>
  void addToHessian(hessian_type & hess,
		    const window_size_t & n,
		    const window_size_t & sbar,
		    const ::pressio::matrixUpperTriangular & hessianTag) const
  {
    for (window_size_t i=0; i < sbar; i++)
    {
      for (window_size_t j=0; j <= i; j++)
      {
        constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();

        auto hess_block =
	  ::pressio::containers::subspan(hess,
					 std::make_pair( (n-i)*romSize_, (n-i+1)*romSize_ ),
					 std::make_pair( (n-j)*romSize_, (n-j+1)*romSize_ ) );

        ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(), one,
				jacobians_.localJacobian( n, i ),
				jacobians_.localJacobian( n, j ),
				one, hess_block);
      }
    }
  }

  template<typename hessian_type>
  void addToHessian(hessian_type & hess,
		    const window_size_t & n,
		    const window_size_t & sbar,
		    const ::pressio::matrixLowerTriangular & hessianTag) const
  {
    for (window_size_t i=0; i < sbar; i++)
    {
      for (window_size_t j=0; j <= i; j++)
      {
        constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
        auto hess_block =
	  ::pressio::containers::subspan(hess,
					 std::make_pair( (n-j)*romSize_, (n-j+1)*romSize_ ),
					 std::make_pair( (n-i)*romSize_, (n-i+1)*romSize_ ) );

        ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(), one,
				jacobians_.localJacobian( n , j ),
				jacobians_.localJacobian( n , i ), one, hess_block);
      }
    }
  }

  // reconstructs fomStateCurrent_ from the stepNum entry of wlsState
  template <typename wls_state_type, typename fom_state_reconstr_t>
  void setCurrentFomState(const wls_state_type & wlsState,
                          const window_size_t & stepNum,
                          const fom_state_reconstr_t & fomStateReconstrObj) const
  {
    const auto wlsCurrentState = ::pressio::containers::span(wlsState,
							     stepNum*romSize_,
							     romSize_);
    fomStateReconstrObj(wlsCurrentState, fomStateCurrent_);
  }

  template <typename time_scheme_t>
  void computeJacobiansOverStencil(const time_scheme_t & timeSchemeObj_,
				   const window_size_t & stepNumGlobal,
           const window_size_t & stepNumLocal,
				   const scalar_t & t,
				   const scalar_t & dt) const
  {

    for (window_size_t i = 0; i < jacStencilSize_; i++)
    {
      auto & jacLocal = jacobians_.localJacobian(stepNumLocal , i );
      timeSchemeObj_.time_discrete_jacobian(fomSystemObj_.fomCRef(),
					    fomStateCurrent_,
					    jacLocal, phi_.get(), t, dt,
					    stepNumGlobal, i);
      Preconditioner(fomSystemObj_.fomCRef(), fomStateCurrent_, jacLocal, t);
    }
  }

  // updates residual and Jacobians at a single time step
  template <
    typename ode_obj_t,
    typename wls_state_type,
    typename fom_state_reconstr_t,
    typename gradient_type>
  void updateResidualAndJacobian(const ode_obj_t & timeSchemeObj_,
				 const wls_state_type & wlsState,
				 const fom_state_reconstr_t & fomStateReconstrObj,
				 const window_size_t & n,
				 const window_size_t & step_s,
				 const scalar_t & ts,
				 scalar_t & rnorm,
				 gradient_type & gradient,
				 const scalar_t & dt) const
  {

    // === reconstruct FOM states ========
    timeSchemeObj_.updateStatesNStep(fomStateCurrent_);
    setCurrentFomState(wlsState,n,fomStateReconstrObj);

    // == Evaluate residual ============
    scalar_t t;
    t = ts + n*dt;
    window_size_t step;
    step = step_s + n;

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("residual");
#endif
    timeSchemeObj_.time_discrete_residual(fomSystemObj_.fomCRef(),
					  fomStateCurrent_,
					  residual_, t, dt, step);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("residual");
#endif

    Preconditioner(fomSystemObj_.fomCRef(), fomStateCurrent_, residual_, t);

    rNormHelper_ = ::pressio::ops::norm2(residual_);
    rnorm += rNormHelper_*rNormHelper_;
    if (innerLoopCounter_ % jacobianUpdateFrequency_ == 0){
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("jacobian");
#endif
      computeJacobiansOverStencil(timeSchemeObj_,step,n,t,dt);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("jacobian");
#endif
    }
  }

private:
  rom_size_t romSize_ = {};
  window_size_t numStepsInWindow_ = {};
  window_size_t jacStencilSize_ = {};
  window_size_t jacobianUpdateFrequency_ = {};

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  ::pressio::rom::impl::FomObjMixin<fom_system_type, true> fomSystemObj_;
#else
  ::pressio::rom::impl::FomObjMixin<fom_system_type, false> fomSystemObj_;
#endif

  std::reference_wrapper<const decoder_jac_t> phi_;

  mutable fom_state_t fomStateCurrent_;
  mutable residual_t residual_;
  const decoder_jac_t J_;

  mutable jacobians_container_t jacobians_;
  const preconditioner_t Preconditioner = {};
  const time_stencil_t timeSchemeObj_;

  mutable window_size_t innerLoopCounter_ = 0;
  mutable scalar_t rNormHelper_ = {};
};

}}}} //end namespace pressio::rom::wls::impl

#endif  // ROM_WLS_IMPL_POLICIES_ROM_WLS_HESSIAN_AND_GRADIENT_SEQUENTIAL_POLICY_IMPL_HPP_
