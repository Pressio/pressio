
#ifndef ROM_GALERKIN_IMPLICIT_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_IMPLICIT_RESIDUAL_POLICY_HPP_

#include "rom_ConfigDefs.hpp"
#include "rom_incremental_solution_base.hpp"
#include "policies/base/ode_implicit_residual_policy_base.hpp"
#include "ode_residual_impl.hpp"

namespace rom{
namespace exp{

template<typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type,
	 typename basis_t>
class romGalerkinImplicitResidualPolicy
  : public ode::policy::implicitResidualPolicyBase<
  romGalerkinImplicitResidualPolicy<state_type, residual_type,
				    model_type, time_type,
				    sizer_type, basis_t>, 1, 0>,
    protected incrementalSolutionBase<romGalerkinImplicitResidualPolicy,
				      state_type, residual_type, model_type,
				      time_type, sizer_type, basis_t>
{

private:
  using base_t = ode::policy::implicitResidualPolicyBase<
  romGalerkinImplicitResidualPolicy<state_type, residual_type,
				    model_type, time_type,
				    sizer_type, basis_t>, 1, 0>;

  using incr_base_t = incrementalSolutionBase<
    romGalerkinImplicitResidualPolicy, state_type, residual_type,
    model_type, time_type, sizer_type, basis_t>;

private:
  using incr_base_t::y0ptr_;
  using incr_base_t::yFull_;
  
public:
  romGalerkinImplicitResidualPolicy(const state_type & y0,
				    const basis_t & phi)
    : incr_base_t(y0), phiPtr_(&phi){}
  
  ~romGalerkinImplicitResidualPolicy() = default;  

private:
  // I need here this auxiliary RHS container 
  // that is needed below to evaluate the spatial residual
  // because the R that we are passed below inthe compute 
  // is supposed to have same size as the 
  // state integrated in time, so for galerkin projection, 
  // we have a modified "space residual" that is 

  //   fnew = \phi^T f(y,...)

  // where f(y,...) is the actual space RHS computed 
  // from the application and so f has to have size suitabel 
  // for the application, whereas fnew has same size of y 

  basis_t const * phiPtr_;
  residual_type appRHS_;
  
private:
  // enable if using types from core package
  template <typename U = state_type,
	    typename T = residual_type,
	    typename std::enable_if<
	      core::meta::is_coreVector<U>::value==true &&
	      core::meta::is_coreVector<T>::value==true
	    >::type * = nullptr>
  void computeImpl(const U & y,
		   T & R,
		   const std::array<U, 1> & oldYs,
		   model_type & model,
		   time_type t,
		   time_type dt)
  {
    // // reconstruct the solution
    // yFull_ = (y0ptr_ + y);
    // *yFull_.data() = phiPtr_->data() * *yFull_.data();
    // std::cout << yFull_.size() << std::endl;
    
    // // eval RHS from target model
    // model.residual(*yFull_.data(), *R.data(), t);

    // ode::impl::implicit_euler_time_discrete_residual(y, oldYs[0], R, dt);
									   

    // // do projection of the space residual
    // // R = phi^T appRHS_
  }  

private:
  friend base_t;
  friend incr_base_t;
  
};//end class

  
}//end namespace exp
}//end namespace rom
#endif 
