
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
	 typename phi_op_t,
	 typename wei_op_t>
class romGalerkinImplicitResidualPolicy
  : public ode::policy::implicitResidualPolicyBase<
  romGalerkinImplicitResidualPolicy<state_type, residual_type,
				    model_type, time_type, sizer_type,
				    phi_op_t, wei_op_t>, 1, 0>
/*protected incrementalSolutionBase<romGalerkinImplicitResidualPolicy,
	      state_type, residual_type, model_type,
	      time_type, sizer_type, phi_op_t, wei_op_t, weighting_op_t>*/
{

private:
  using base_t = ode::policy::implicitResidualPolicyBase<
  romGalerkinImplicitResidualPolicy<state_type, residual_type,
				    model_type, time_type,
				    sizer_type, phi_op_t, wei_op_t>, 1, 0>;

  // using incr_base_t = incrementalSolutionBase<
  //   romGalerkinImplicitResidualPolicy, state_type, residual_type,
  //   model_type, time_type, sizer_type, phi_op_t, wei_op_t>;
// private:
//   using incr_base_t::y0ptr_;
//   using incr_base_t::yFull_;
  
private:
  phi_op_t * phiOp_;
  wei_op_t * WOp_;
  residual_type appRHS_;
  state_type yFOM_;
  
public:
  romGalerkinImplicitResidualPolicy(const state_type & y0fom,
				    const state_type & y0r,
				    phi_op_t & phiOp,
				    wei_op_t & WOp)
    : phiOp_(&phiOp), WOp_(&WOp), appRHS_(y0fom.size()), yFOM_(y0fom.size())
  {}
  
  ~romGalerkinImplicitResidualPolicy() = default;  

private:
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
    // reconstruct the reduced actual solution
    //    yFull_ = y;//(*y0ptr_ + y);
    // std::cout << *yFull_.data() << std::endl;
    // std::cout << "--------------" << std::endl;
    // std::cout << *y0ptr_->data() << std::endl;
    // std::cout << "--------------" << std::endl;
    // std::cout << *y.data() << std::endl;

    // reconstruct FOM state vector for computing space residual
    phiOp_->leftMultiply(y, yFOM_);

    // eval RHS from target model
    model.residual(*yFOM_.data(), *appRHS_.data(), t);
    
    if (sizer_type::getSize(R)==0)
      sizer_type::matchSize(y, R);

    WOp_->leftMultiply(appRHS_, R);

    // *yFOM_.data() = *phiPtr_->data() * (*y.data());    
    // // eval RHS from target model
    // model.residual(*yFOM_.data(), *appRHS_.data(), t);
    // if (R.empty())
    //   R.resize(y.size());

    // *R.data() = *phiTPtr_->data() * (*appRHS_.data());
    ode::impl::implicit_euler_time_discrete_residual(y, oldYs[0], R, dt);
  }  

private:
  friend base_t;
  //friend incr_base_t;
  
};//end class

  
}//end namespace exp
}//end namespace rom
#endif 
