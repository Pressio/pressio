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

#ifndef ROM_WLS_HESSIAN_GRADIENT_SEQUENTIAL_POLICY_IMPL_HPP_
#define ROM_WLS_HESSIAN_GRADIENT_SEQUENTIAL_POLICY_IMPL_HPP_

#include "rom_wls_preconditioners_impl.hpp" // this needs to get taken out of the impl namespace
#include "rom_wls_jacobians_container.hpp"

/*
This header file contains the class used for computing the hessian and gradients in the WLS system.
The hessian_gradient policy is responsible for assembling Jw^T Jw and J^T rw,
where Jw is the reduced windowed Jacobian and rw the residual.
The policy initializes in memory:
wlsJacs: this is a vector container of n_s -1 local Jacobians, J, where n_s is the width of the time stencil (e.g., n_s = 3 for BDF2)
residual: this is the residual vector for the FOM
yFOM_current: this is a working variable for the fom state.
*/

namespace pressio{ namespace rom{ namespace wls{ namespace impl{

template<
  typename fom_type,
  typename decoder_t,
  typename hessian_matrix_structure_tag,
  typename preconditioner_t = ::pressio::rom::wls::preconditioners::impl::NoPreconditioner,
  typename jacobians_container_t = typename ::pressio::rom::wls::impl::nonFrozenJacobiansContainer<typename decoder_t::jacobian_type>
  >
class HessianGradientSequentialPolicy
{
public:
  using scalar_t                = typename fom_type::scalar_type;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using decoder_jac_t           = typename decoder_t::jacobian_type;
  using residual_t              = fom_state_t;
  /* wls_jacs_t is for now a vector of item type = decoder_jac_t */
  using wls_jacs_t	      = std::vector<decoder_jac_t>;
public:
  HessianGradientSequentialPolicy(rom_size_t romSize,
				  window_size_t numStepsInWindow,
				  const decoder_t & decoderObj,
				  const fom_type & appObj,
          const fom_state_t & fomState,
          window_size_t timeStencilSize,
          const window_size_t jacobianUpdateFrequency = 1)
    : romSize_(romSize),
      phi_(decoderObj.getReferenceToJacobian()),
      appObj_(appObj),
      fomStateCurrent_(fomState),
      jacStencilSize_(std::min(timeStencilSize+1, numStepsInWindow)),
      // construct wls Jacobians from jacobian of the decoder: we might need to change this later
      Jacobians_( timeStencilSize, numStepsInWindow , decoderObj.getReferenceToJacobian()), 
      residual_( appObj.velocity( *fomState.data(), ::pressio::utils::constants::zero<scalar_t>()) ),
      jacobianUpdateFrequency_(jacobianUpdateFrequency)
  {
    if(  (std::is_same<jacobians_container_t,::pressio::rom::wls::impl::nonFrozenJacobiansContainer<typename decoder_t::jacobian_type> >::value) and 
         ( jacobianUpdateFrequency != 1 ) ){
      ::pressio::utils::io::print_stdout("Error, using nonFrozenJacobianContainer with jacobianUpdateFrequency > 1. To run with jacobianUpdateFrequency > 1, use frozenJacobiansContainer  \n");
      ::pressio::utils::io::print_stdout("Setting jacobianUpdateFrequency = 1  \n");
      jacobianUpdateFrequency_ = 1;
    }
  }

public:
  template <
    typename ode_obj_t,
    typename wls_state_type,
    typename fom_state_reconstr_t,
    typename hess_type,
    typename gradient_type
  >
  void operator()(ode_obj_t & timeSchemeObj,
                  const wls_state_type  & wlsState,
                  const wls_state_type & wlsStateIC,
                  hess_type & hess,
                  gradient_type & gradient,
                  const fom_state_reconstr_t & fomStateReconstrObj,
                  const scalar_t dt,
                  const window_size_t numStepsInWindow,
                  const scalar_t ts ,
                  const window_size_t step_s,
                  scalar_t & rnorm) const
  {
    constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
    constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();

    window_size_t stepNumLocal = 0;
    scalar_t t = ts + stepNumLocal*dt;
    window_size_t stepNumGlobal = step_s + stepNumLocal;
    //if (innerLoopCounter_%freezeFreq_ == 0){
    ::pressio::ops::set_zero(hess);
    //}
    ::pressio::ops::set_zero(gradient);

    //get access to the state at the first window
    setCurrentFomState(wlsState, 0, fomStateReconstrObj);

    //reconstruct the FOM states from the previous window/ICs
    timeSchemeObj.updateStatesFirstStep(wlsStateIC, fomStateReconstrObj);

    //compute the time discrete residual
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("residual");
#endif
    timeSchemeObj.time_discrete_residual(appObj_, fomStateCurrent_, residual_, ts, dt, stepNumGlobal);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("residual");
#endif
    Preconditioner(appObj_, fomStateCurrent_, residual_, t);

    //increment the norm
    rnorm += ::pressio::ops::norm2(residual_);

    // compute jacobian over stencil
    if (innerLoopCounter_%jacobianUpdateFrequency_ == 0){
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("jacobian");
#endif
      computeJacobiansOverStencil(timeSchemeObj, stepNumGlobal, stepNumLocal , t,dt);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("jacobian");
#endif
    } //end if for freezing

    // add to local block of hessian
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("hessian computation");
#endif
    auto hess_block = ::pressio::containers::subspan(hess,
						     std::make_pair( stepNumLocal*romSize_,(stepNumLocal+1)*romSize_ ) ,
						     std::make_pair( stepNumLocal*romSize_,(stepNumLocal+1)*romSize_) );
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			    one, Jacobians_.getLocalJacobian(stepNumLocal,0),
			    Jacobians_.getLocalJacobian(stepNumLocal,0), zero, hess_block);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("hessian computation");
#endif

    // compute gradient[n*romSize_:(n+1)*romSize] += J^T r
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("gradient computation");
#endif
    auto gradientView = ::pressio::containers::span(gradient, stepNumLocal*romSize_, romSize_);
    ::pressio::ops::product(::pressio::transpose(),
			    one, Jacobians_.getLocalJacobian(stepNumLocal,0), residual_,
			    one, gradientView);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("gradient computation");
#endif
    for (window_size_t stepNumLocal = 1; stepNumLocal < numStepsInWindow; stepNumLocal++)
    {
      updateResidualAndJacobian(timeSchemeObj, wlsState, fomStateReconstrObj, stepNumLocal,
				step_s, ts, rnorm, gradient, dt);

      const window_size_t sbar = std::min(stepNumLocal, jacStencilSize_);
      for (window_size_t i=0; i < sbar; i++)
      {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
        timer->start("gradient computation");
#endif
      	auto gradientView = ::pressio::containers::span(gradient, (stepNumLocal-i)*romSize_, romSize_);
        ::pressio::ops::product(::pressio::transpose(),
				one, Jacobians_.getLocalJacobian( stepNumLocal, i ), residual_,
				one, gradientView);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
        timer->stop("gradient computation");
#endif
      }

//      if (innerLoopCounter_%freezeFreq_ == 0){  
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("hessian computation");
#endif
      addToHessian(hess, stepNumLocal, sbar, hessian_matrix_structure_tag());
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("hessian computation");
#endif
//      }

    }//end loop over stepsInWindow
  innerLoopCounter_ += 1;
  }//end operator()

private:
  template<typename hessian_type>
  void addToHessian(hessian_type & hess,
		    const window_size_t & n,
		    const window_size_t & sbar,
		    const ::pressio::matrixUpperTriangular & hessianTag) const
  {
    for (window_size_t i=0; i < sbar; i++){
      for (window_size_t j=0; j <= i; j++){
        constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
        auto hess_block = ::pressio::containers::subspan(hess,
							 std::make_pair( (n-i)*romSize_, (n-i+1)*romSize_ ),
							 std::make_pair( (n-j)*romSize_,(n-j+1)*romSize_ ) );
        ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(), one,
				Jacobians_.getLocalJacobian( n, i ),
				Jacobians_.getLocalJacobian( n, j ), one, hess_block);
      }
    }// end assembling local component of global Hessian
  }

  template<typename hessian_type>
  void addToHessian(hessian_type & hess,
		    const window_size_t & n,
		    const window_size_t & sbar,
		    const ::pressio::matrixLowerTriangular & hessianTag) const
  {
    for (window_size_t i=0; i < sbar; i++){
      for (window_size_t j=0; j <= i; j++)
      {
        constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
        auto hess_block = ::pressio::containers::subspan(hess,
							 std::make_pair( (n-j)*romSize_, (n-j+1)*romSize_ ),
							 std::make_pair( (n-i)*romSize_,(n-i+1)*romSize_ ) );

        ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
				one,
				Jacobians_.getLocalJacobian( n , j ),
				Jacobians_.getLocalJacobian( n , i ), one, hess_block);
      }
    }// end assembling local component of global Hessian
  }

  // reconstructs fomStateCurrent_ from the stepNum entry of wlsState
  template <typename wls_state_type, typename fom_state_reconstr_t>
  void setCurrentFomState(const wls_state_type & wlsState,
                          const window_size_t & stepNum,
                          const fom_state_reconstr_t & fomStateReconstrObj) const
  {
    const auto wlsCurrentState = ::pressio::containers::span(wlsState, stepNum*romSize_, romSize_);
    fomStateReconstrObj(wlsCurrentState, fomStateCurrent_);
  }

  template <typename time_scheme_t>
  void computeJacobiansOverStencil(const time_scheme_t & timeSchemeObj,
				   const window_size_t & stepNumGlobal,
           const window_size_t & stepNumLocal,
				   const scalar_t & t,
				   const scalar_t & dt) const
  {
    for (window_size_t i = 0; i < jacStencilSize_; i++)
    {
      auto & jacLocal = Jacobians_.getLocalJacobian(stepNumLocal , i );
      timeSchemeObj.time_discrete_jacobian(appObj_,
                                           fomStateCurrent_,
                                           jacLocal,
                                           phi_, t, dt, stepNumGlobal, i);
      //if (timeSchemeObj.jacobianNeedsRecomputing(jacStencilSize_ -i -1)){
      Preconditioner(appObj_, fomStateCurrent_, jacLocal, t);
      //}
    }
  }

  // updates residual and Jacobians at a single time step
  template <
    typename ode_obj_t,
    typename wls_state_type,
    typename fom_state_reconstr_t,
    typename gradient_type>
  void updateResidualAndJacobian(const ode_obj_t & timeSchemeObj,
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
      timeSchemeObj.updateStatesNStep(fomStateCurrent_);
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
      timeSchemeObj.time_discrete_residual(appObj_, fomStateCurrent_, residual_, t, dt, step);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("residual");
#endif

      Preconditioner(appObj_,fomStateCurrent_,residual_,t);

      rnorm += ::pressio::ops::norm2(residual_);

      if (innerLoopCounter_%jacobianUpdateFrequency_ == 0){
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("jacobian");
#endif
      computeJacobiansOverStencil(timeSchemeObj,step,n,t,dt);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("jacobian");
#endif
      }
  }

private:
  rom_size_t romSize_;
  const decoder_jac_t & phi_;
  const fom_type & appObj_;
  mutable fom_state_t fomStateCurrent_; //working variable for the FOM state
  window_size_t jacStencilSize_;

  mutable residual_t residual_;	        // working variable for the time discrete residual
  const preconditioner_t Preconditioner{};
  window_size_t jacobianUpdateFrequency_;
  mutable window_size_t innerLoopCounter_ = 0;
  jacobians_container_t Jacobians_;

};

}}}} //end namespace pressio::rom::wls::impl
#endif





// template <
//   typename ode_obj_t,
//   typename wls_state_type,
//   typename fom_state_reconstr_t,
//   typename aux_states_container_t>
// void computeResidualNorm(const fom_type & appObj,
// 			   ode_obj_t & odeObj_,
// 			   const wls_state_type  & wlsState,
// 			   const wls_state_type & wlsStateIC,
// 			   const fom_state_reconstr_t & fomStateReconstrObj_,
// 			   const scalar_t dt,
// 			   const std::size_t numStepsInWindow,
// 			   const scalar_t ts ,
// 			   aux_states_container_t & auxStatesContainer,
// 			   const int step_s,
// 			   scalar_t & rnorm) const
// {
//   rnorm = ::pressio::utils::constants::zero<scalar_t>();
//   int n = 0;
//   auto t = ts + n*dt;
//   int step = step_s + n;

//   const auto wlsCurrentState = ::pressio::containers::span(wlsState,0,romSize_);
//   fomStateReconstrObj_(wlsCurrentState,fomStateCurrent_);
//   odeObj_.updateStatesFirstStep(wlsStateIC,fomStateReconstrObj_);

//   odeObj_.time_discrete_residual(appObj,fomStateCurrent_,residual_,ts,dt,step);
//   rnorm += ::pressio::ops::norm2(residual_);

//   for (int n = 1; n < numStepsInWindow; n++)
//   {
//     // === reconstruct FOM states ========
//     odeObj_.updateStatesNStep(fomStateCurrent_);
//     const auto wlsCurrentState = ::pressio::containers::span(wlsState,n*romSize_,romSize_);
//     fomStateReconstrObj_(wlsCurrentState,fomStateCurrent_);

//     // == Evaluate residual ============
//     t = ts + n*dt;
//     step = step_s + n;
//     odeObj_.time_discrete_residual(appObj,fomStateCurrent_,residual_,t,dt,step);
//     rnorm += ::pressio::ops::norm2(residual_);
//   }
// }//end computeResidualNorm
