
#ifndef ROM_WLS_HESSIAN_GRADIENT_SEQUENTIAL_POLICY_IMPL_HPP_
#define ROM_WLS_HESSIAN_GRADIENT_SEQUENTIAL_POLICY_IMPL_HPP_

/*
This header file contains the class object used for computing the hessian and gradients in the WLS system.
The hessian_gradient policy is responsible for assembling Jw^T Jw and J^T rw,
where Jw is the reduced windowed Jacobian and rw the residual.
The policy initializes in memory:
wlsJacs: this is a vector container of n_s -1 local Jacobians, J, where n_s is the width of the time stencil (e.g., n_s = 3 for BDF2)
residual: this is the residual vector for the FOM
yFOM_current: this is a working variable for the fom state.
*/

namespace pressio{ namespace rom{ namespace wls{ namespace impl{

// Helper function to do c += A^T b
template <typename Mat, typename Vec>
void local_vec_update(const Mat & A,
		      const  Vec & b,
		      Vec  & cOut,
		      const int & scol,
		      const int & colSize)
{
  auto gradientView = ::pressio::containers::span(cOut,scol,colSize);
  auto tmp = pressio::containers::ops::dot(A, b);
  for (int k =0; k< colSize; k++){
    gradientView[k] += tmp[k];}
}



template<typename fom_type, typename decoder_t>
class HessianGradientSequentialPolicy
{

private:
  using scalar_t                = typename fom_type::scalar_type;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using decoder_jac_t           = typename decoder_t::jacobian_t;

  // currently have this as a vector of jacobians, can change to container
  /* For now, set wls_jacs_t to be a vector of item type = decoder_jac_t
   * this is a similar assumption as done for lspg
   */
  using wls_jacs_t		= std::vector<decoder_jac_t>;

  //maybe we need to change this for hyper reduction?
  using residual_t              = fom_state_t;

  mutable wls_jacs_t wlsJacs_;
  mutable fom_state_t yFOM_current_; //working variable for the FOM state
  mutable residual_t residual_;	     // working variable for the time discrete residual
  int romSize_;
  int fomSize_;
  int time_stencil_size;
  const decoder_jac_t & phi_;

public:
  HessianGradientSequentialPolicy(const fom_type & appObj,
				  const fom_state_t & yFOM,
				  int numStepsInWindow,
				  int time_stencil_size,
				  const decoder_t & decoderObj)
    : // construct wls Jacobians from jacobian of the decoder: we might need to change this later
      wlsJacs_(time_stencil_size,decoderObj.getReferenceToJacobian()),
      phi_(decoderObj.getReferenceToJacobian()),
      residual_( appObj.velocity( *yFOM.data() , ::pressio::utils::constants::zero<scalar_t>()) ),
      yFOM_current_(yFOM)
  {
    this->romSize_ = decoderObj.getReferenceToJacobian().numVectors();
    this->time_stencil_size = time_stencil_size;
    this->fomSize_ = yFOM.size();
  }


  template <
    typename ode_obj_t,
    typename wls_state_type,
    typename fom_state_reconstr_t,
    typename hess_type,
    typename gradient_type>
  void operator()(const fom_type & appObj,
                  ode_obj_t & odeObj_,
                  const wls_state_type  & wlsState,
                  const wls_state_type & wlsStateIC,
                  hess_type & hess,
                  gradient_type & gradient,
                  const fom_state_reconstr_t & fomStateReconstrObj_,
                  const scalar_t dt,
                  const std::size_t numStepsInWindow,
                  const scalar_t ts ,
                  const int step_s,
                  scalar_t & rnorm) const
  {
    int n = 0;
    scalar_t t = ts + n*dt;
    int step = step_s + n;
    hess.setZero(); //ensure zero as we use +=
    gradient.setZero();

    //get access to the state at the first window
    const auto wlsCurrentState = ::pressio::containers::span(wlsState, 0, romSize_);
    fomStateReconstrObj_(wlsCurrentState, yFOM_current_);
    //reconstruct the FOM states from the previous window/ICs
    odeObj_.updateStatesFirstStep(wlsStateIC, fomStateReconstrObj_);

    //compute the time discrete residual
    odeObj_.time_discrete_residual(appObj, yFOM_current_, residual_, ts, dt, step);
    //increment the norm
    rnorm += ::pressio::containers::ops::norm2(residual_);

    // this is a temporary work around for the += issue.
    hess_type C(romSize_,romSize_);
    for (int i = 0; i < time_stencil_size; i++){
      //compute the Jacobian for each state in the stencil
      odeObj_.time_discrete_jacobian(appObj, yFOM_current_,
				     wlsJacs_[time_stencil_size - i - 1],
				     phi_, n*dt, dt, step, i);
    }
    // add to local block of hessian
    auto hess_block = ::pressio::containers::subspan( hess,
						      std::make_pair( n*romSize_,(n+1)*romSize_ ) ,
						      std::make_pair( n*romSize_,(n+1)*romSize_) );
    ::pressio::containers::ops::dot(wlsJacs_[time_stencil_size-1], wlsJacs_[time_stencil_size-1], hess_block);

    // compute gradient[n*romSize_:(n+1)*romSize] += J^T r
    local_vec_update(wlsJacs_[time_stencil_size-1], residual_, gradient, n*romSize_, romSize_);

    for (int n = 1; n < numStepsInWindow; n++)
    {
      // === reconstruct FOM states ========
      odeObj_.updateStatesNStep(yFOM_current_);
      const auto wlsCurrentState = ::pressio::containers::span(wlsState, n*romSize_, romSize_);
      fomStateReconstrObj_(wlsCurrentState, yFOM_current_);

      // == Evaluate residual ============
      t = ts + n*dt;
      step = step_s + n;
      odeObj_.time_discrete_residual(appObj, yFOM_current_, residual_, t, dt, step);
      rnorm += ::pressio::containers::ops::norm2(residual_);

      for (int i = 0; i < time_stencil_size; i++){
	//compute the Jacobian for each state in the stencil
	odeObj_.time_discrete_jacobian(appObj, yFOM_current_,
				       wlsJacs_[time_stencil_size - i - 1],
				       phi_, n*dt, dt, step, i);
      }

      // == Update everything
      int sbar = std::min(n,time_stencil_size);
      for (int i=0; i < sbar; i++)
      {
	local_vec_update(wlsJacs_[time_stencil_size-i-1], residual_, gradient, (n-i)*romSize_, romSize_);
      }

      // == Assemble local component of global Hessian //
      for (int i=0; i < sbar; i++){
	for (int j=0; j <= i; j++)
	{
	  auto hess_block = ::pressio::containers::subspan(hess,
							   std::make_pair( (n-i)*romSize_, (n-i+1)*romSize_ ),
							   std::make_pair( (n-j)*romSize_,(n-j+1)*romSize_ ) );
	  ::pressio::containers::ops::dot(wlsJacs_[time_stencil_size-i-1], wlsJacs_[time_stencil_size-j-1], C);

	  for (int k = 0; k<romSize_; k++){
	    for (int l = 0; l<romSize_; l++){
	      hess_block(k,l) = hess_block(k,l) + C(k,l);
	    }
	  }

	  auto hess_block2 = ::pressio::containers::subspan(hess,
							    std::make_pair( (n-j)*romSize_, (n-j+1)*romSize_),
							    std::make_pair( (n-i)*romSize_,(n-i+1)*romSize_) );
	  for (int k = 0; k<romSize_; k++){
	    for (int l = 0; l<romSize_; l++){
	      hess_block2(l,k) = hess_block(k,l);
	    }
	  }
	}
      }// end assembling local component of global Hessian
    }//end loop over stepsInWindow

  }//end operator()


  // FRizzi: is this used? I commented it out and things work anyway
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
  //   fomStateReconstrObj_(wlsCurrentState,yFOM_current_);
  //   odeObj_.updateStatesFirstStep(wlsStateIC,fomStateReconstrObj_);

  //   odeObj_.time_discrete_residual(appObj,yFOM_current_,residual_,ts,dt,step);
  //   rnorm += ::pressio::containers::ops::norm2(residual_);

  //   for (int n = 1; n < numStepsInWindow; n++)
  //   {
  //     // === reconstruct FOM states ========
  //     odeObj_.updateStatesNStep(yFOM_current_);
  //     const auto wlsCurrentState = ::pressio::containers::span(wlsState,n*romSize_,romSize_);
  //     fomStateReconstrObj_(wlsCurrentState,yFOM_current_);

  //     // == Evaluate residual ============
  //     t = ts + n*dt;
  //     step = step_s + n;
  //     odeObj_.time_discrete_residual(appObj,yFOM_current_,residual_,t,dt,step);
  //     rnorm += ::pressio::containers::ops::norm2(residual_);
  //   }
  // }//end computeResidualNorm

};

}}}}
#endif
