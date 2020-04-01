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

#include "rom_wls_preconditioners_impl.hpp"

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
  typename preconditioner_t = ::pressio::rom::wls::preconditioners::impl::NoPreconditioner
  >
class HessianGradientSequentialPolicy
{
  using scalar_t                = typename fom_type::scalar_type;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using decoder_jac_t           = typename decoder_t::jacobian_type;
  using residual_t              = fom_state_t;

  // currently have this as a vector of jacobians, can change to container
  /* For now, set wls_jacs_t to be a vector of item type = decoder_jac_t
   * this is a similar assumption as done for lspg
   */
  using wls_jacs_t	      = std::vector<decoder_jac_t>;

public:
  HessianGradientSequentialPolicy(const fom_type & appObj,
                                  const fom_state_t & fomState,
                                  window_size_t numStepsInWindow,
                                  window_size_t time_stencil_size,
                                  const decoder_t & decoderObj,
                                  rom_size_t romSize)
    : // construct wls Jacobians from jacobian of the decoder: we might need to change this later
      appObj_(appObj),
      jac_stencil_size_(std::min(time_stencil_size+1, numStepsInWindow)),
      wlsJacs_(std::min(time_stencil_size+1,numStepsInWindow), decoderObj.getReferenceToJacobian()),
      phi_(decoderObj.getReferenceToJacobian()),
      residual_( appObj.velocity( *fomState.data() , ::pressio::utils::constants::zero<scalar_t>()) ),
      fomStateCurrent_(fomState)
  {
    this->romSize_ = romSize;
    this->time_stencil_size_ = time_stencil_size;
  }


public:
  template <
    typename ode_obj_t,
    typename wls_state_type,
    typename fom_state_reconstr_t,
    typename hess_type,
    typename gradient_type>
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

    window_size_t n = 0;
    scalar_t t = ts + n*dt;
    window_size_t step = step_s + n;
    ::pressio::ops::set_zero(hess);
    ::pressio::ops::set_zero(gradient);

    //get access to the state at the first window
    setCurrentFomState(wlsState,0,fomStateReconstrObj);

    //reconstruct the FOM states from the previous window/ICs
    timeSchemeObj.updateStatesFirstStep(wlsStateIC, fomStateReconstrObj);

    //compute the time discrete residual
    timeSchemeObj.time_discrete_residual(appObj_, fomStateCurrent_, residual_, ts, dt, step);
    Preconditioner(appObj_,fomStateCurrent_,residual_,t);

    //increment the norm
    rnorm += ::pressio::ops::norm2(residual_);

    // compute jacobian over stencil
    computeJacobiansOverStencil(timeSchemeObj, step,t,dt);

    // add to local block of hessian
    auto hess_block = ::pressio::containers::subspan( hess,
						      std::make_pair( n*romSize_,(n+1)*romSize_ ) ,
						      std::make_pair( n*romSize_,(n+1)*romSize_) );
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
					one, wlsJacs_[jac_stencil_size_-1],
					wlsJacs_[jac_stencil_size_-1], zero, hess_block);

    // compute gradient[n*romSize_:(n+1)*romSize] += J^T r
    auto gradientView = ::pressio::containers::span(gradient, n*romSize_, romSize_);
    ::pressio::ops::product(::pressio::transpose(),
					one, wlsJacs_[jac_stencil_size_-1], residual_,
					one, gradientView);

    for (window_size_t n = 1; n < numStepsInWindow; n++)
    {
      updateResidualAndJacobian(timeSchemeObj, wlsState, fomStateReconstrObj, n,
				step_s, ts, rnorm, gradient, dt);

      const window_size_t sbar = std::min(n, jac_stencil_size_);
      for (window_size_t i=0; i < sbar; i++)
      {
      	auto gradientView = ::pressio::containers::span(gradient, (n-i)*romSize_, romSize_);
        ::pressio::ops::product(::pressio::transpose(),
					    one, wlsJacs_[jac_stencil_size_-i-1], residual_,
					    one, gradientView);
      }
      addToHessian(hess, n, sbar, hessian_matrix_structure_tag());
    }//end loop over stepsInWindow
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
        ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
				one,
				wlsJacs_[jac_stencil_size_-i-1],
				wlsJacs_[jac_stencil_size_-j-1], one, hess_block);
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
				wlsJacs_[jac_stencil_size_-j-1],
				wlsJacs_[jac_stencil_size_-i-1], one, hess_block);
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
				   const window_size_t & stepNum,
				   const scalar_t & t,
				   const scalar_t & dt) const
  {
    for (window_size_t i = 0; i < jac_stencil_size_; i++)
    {
      timeSchemeObj.time_discrete_jacobian(appObj_,
                                           fomStateCurrent_,
                                           wlsJacs_[jac_stencil_size_ -i -1],
                                           phi_, t, dt, stepNum, i);

      if (timeSchemeObj.jacobianNeedsRecomputing(jac_stencil_size_ -i -1)){
        Preconditioner(appObj_, fomStateCurrent_, wlsJacs_[jac_stencil_size_ -i -1], t);
      }
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
      timeSchemeObj.time_discrete_residual(appObj_, fomStateCurrent_, residual_, t, dt, step);
      Preconditioner(appObj_,fomStateCurrent_,residual_,t);
      rnorm += ::pressio::ops::norm2(residual_);
      computeJacobiansOverStencil(timeSchemeObj,step,t,dt);
  }

private:
  mutable wls_jacs_t wlsJacs_;
  mutable fom_state_t fomStateCurrent_; //working variable for the FOM state
  mutable residual_t residual_;	        // working variable for the time discrete residual
  rom_size_t romSize_;
  window_size_t time_stencil_size_;
  window_size_t jac_stencil_size_;
  const decoder_jac_t & phi_;
  const fom_type & appObj_;
  const preconditioner_t Preconditioner{};
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
