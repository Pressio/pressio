
#ifndef ROM_GALERKIN_IMPLICIT_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_IMPLICIT_RESIDUAL_POLICY_HPP_

#include "../rom_ConfigDefs.hpp"
#include "../../../ode/src/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../../../ode/src/ode_residual_impl.hpp"
//#include "rom_incremental_solution_base.hpp"
#include "../../../CORE_ALL"

namespace rompp{
namespace rom{
namespace exp{

template<typename state_type, 
	 typename residual_type,
	 typename model_type,
	 typename phi_type,
	 typename A_type>
class RomGalerkinImplicitResidualPolicy
  : public ode::policy::ImplicitResidualPolicyBase<
  RomGalerkinImplicitResidualPolicy<state_type, residual_type,
				    model_type, phi_type, A_type>, 1, 0>
{

private:
  using base_t = ode::policy::ImplicitResidualPolicyBase<
  RomGalerkinImplicitResidualPolicy<state_type, residual_type,
				    model_type, phi_type, A_type>, 1, 0>;
private:
  state_type yFOM_;
  residual_type appRHS_;
  phi_type * phi_;
  A_type * A_;
  
public:
  RomGalerkinImplicitResidualPolicy(const state_type & y0fom,
				    const residual_type & r0fom,
				    phi_type & phiOp,
				    A_type & AOp)
    : yFOM_(y0fom), appRHS_(r0fom), phi_(&phiOp), A_(&AOp){}
  
  ~RomGalerkinImplicitResidualPolicy() = default;  

private:
  template <typename U = state_type,
	    typename T = residual_type,
	    typename scalar_type,
	    typename std::enable_if<
	      core::meta::is_core_vector_wrapper<U>::value==true &&
	      core::meta::is_core_vector_wrapper<T>::value==true
	    >::type * = nullptr>
  void computeImpl(const U & y,
		   T & R,
		   const std::array<U, 1> & oldYs,
		   model_type & model,
		   scalar_type t,
		   scalar_type dt)
  {
    // y coming in is the REDUCED state, so we need to reconstruct full state
    // assert( y.globalSize() == phi_->globalCols() );
    // yFOM_ = core::matrixVectorProduct(*phi_, y);
    // std::cout << "yFOM_ = " << yFOM_.globalSize() << std::endl;

    // // reconstruct the reduced actual solution
    // //    yFull_ = y;//(*y0ptr_ + y);
    // // std::cout << *yFull_.data() << std::endl;
    // // std::cout << "--------------" << std::endl;
    // // std::cout << *y0ptr_->data() << std::endl;
    // // std::cout << "--------------" << std::endl;
    // // std::cout << *y.data() << std::endl;

    // // reconstruct FOM state vector for computing space residual
    // phiOp_->leftMultiply(y, yFOM_);

    // // eval RHS from target model
    // model.residual(*yFOM_.data(), *appRHS_.data(), t);
    
    // if (sizer_type::getSize(R)==0)
    //   sizer_type::matchSize(y, R);

    // WOp_->leftMultiply(appRHS_, R);

    // // *yFOM_.data() = *phiPtr_->data() * (*y.data());    
    // // // eval RHS from target model
    // // model.residual(*yFOM_.data(), *appRHS_.data(), t);
    // // if (R.empty())
    // //   R.resize(y.size());

    // // *R.data() = *phiTPtr_->data() * (*appRHS_.data());
    // ode::impl::implicit_euler_time_discrete_residual(y, oldYs[0], R, dt);
  }  

private:
  friend base_t;  
};//end class
  
}//end namespace exp
}//end namespace rom
}//end namespace rompp
#endif 
