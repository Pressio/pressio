
#ifndef ROM_GALERKIN_IMPLICIT_JACOBIAN_POLICY_HPP_
#define ROM_GALERKIN_IMPLICIT_JACOBIAN_POLICY_HPP_

#include "rom_ConfigDefs.hpp"
#include "policies/base/ode_jacobian_policy_base.hpp"
#include "ode_jacobian_impl.hpp"
//#include "rom_incremental_solution_base.hpp"

namespace rom{
namespace exp{

template<typename state_type,
	 typename jacobian_type,
	 typename model_type,
	 typename sizer_type,
	 typename phi_type,
	 typename A_type>
class RomGalerkinImplicitJacobianPolicy 
  : public ode::policy::JacobianPolicyBase<
  RomGalerkinImplicitJacobianPolicy<state_type, jacobian_type,
				    model_type, sizer_type,
				    phi_type, A_type> >
{

private:
  using base_t = ode::policy::JacobianPolicyBase<
  RomGalerkinImplicitJacobianPolicy<state_type, jacobian_type,
				    model_type, sizer_type,
				    phi_type, A_type> >;
private:
  state_type yFOM_;
  jacobian_type JFOM_;
  phi_type * phi_;
  A_type * A_;
  
public:
  RomGalerkinImplicitJacobianPolicy(const state_type & y0fom,
				    const jacobian_type & j0fom,
				    phi_type & phiOp,
				    A_type & AOp)
    : yFOM_(y0fom), JFOM_(j0fom), phi_(&phiOp), A_(&AOp){}

  ~RomGalerkinImplicitJacobianPolicy() = default;
  
private:
  template <typename U = state_type,
	    typename T = jacobian_type,
	    typename scalar_type,
	    typename std::enable_if<
	      core::meta::is_core_vector<U>::value==true &&
	      core::meta::is_core_matrix_wrapper<T>::value==true
	    >::type * = nullptr>
  void computeImpl(const U & y,
		   T & J,
		   model_type & model,
		   scalar_type t,
		   scalar_type dt)
  {
    //
  }  

private:
  friend base_t;
};//end class
  
}//end namespace exp
}//end namespace rom
#endif 






// // reconstruct the solution
// //    yFull_ = y;//*y0ptr_ + y;

// if (J.rows()==0)
//   J.resize(y.size(), y.size());

// // reconstruct FOM state vector for computing space residual
// phiOp_->leftMultiply(y, yFOM_);

// // first eval space jacobian
// model.jacobian( *yFOM_.data(), *JFOM_.data(), t);

// WOp_->leftRightMultiply(JFOM_, J);

// // update from time discrete residual
// ode::impl::implicit_euler_time_discrete_jacobian(J, dt);
    
// //////

// *yFOM_.data() = *phiPtr_->data() * (*y.data());
// if (J2_.rows()==0)
//   J2_.resize(yFOM_.size(), yFOM_.size());    
// // first eval space jac
// model.jacobian( *yFOM_.data(), *J2_.data(), t);
// Eigen::MatrixXd JJ(*J2_.data());
// auto a = phiTPtr_->data();
// auto b = phiPtr_->data();
// auto res = (*a) * JJ * (*b);
// // std::cout << "resSize "
// // 	      << a->rows() << " " << a->cols() << " "
// // 	      << JJ.rows() << " " << JJ.cols() << " "
// // 	      << b->rows() << " " << b->cols() << " "
// // 	      << std::endl;
// (*J.data()) = res.sparseView();    
// // update from time discrete residual
// ode::impl::implicit_euler_time_discrete_jacobian(J, dt);
