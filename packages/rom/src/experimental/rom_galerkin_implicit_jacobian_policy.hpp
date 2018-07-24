
#ifndef ROM_GALERKIN_IMPLICIT_JACOBIAN_POLICY_HPP_
#define ROM_GALERKIN_IMPLICIT_JACOBIAN_POLICY_HPP_

#include "rom_ConfigDefs.hpp"
#include "rom_incremental_solution_base.hpp"
#include "policies/base/ode_jacobian_policy_base.hpp"
#include "ode_jacobian_impl.hpp"

namespace rom{
namespace exp{

template<typename state_type,
	 typename jacobian_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type,
	 typename phi_op_t,
	 typename wei_op_t>
class romGalerkinImplicitJacobianPolicy 
  : public ode::policy::JacobianPolicyBase<
  romGalerkinImplicitJacobianPolicy<state_type, jacobian_type,
				    model_type, time_type,
				    sizer_type, phi_op_t, wei_op_t> >
    // protected incrementalSolutionBase<romGalerkinImplicitJacobianPolicy,
    // 				      state_type, jacobian_type,
    // 				      model_type, time_type, sizer_type,
    // 				      phi_op_t, wei_op_t>
{

private:
  using base_t = ode::policy::JacobianPolicyBase<
  romGalerkinImplicitJacobianPolicy<state_type, jacobian_type,
				    model_type, time_type,
				    sizer_type, phi_op_t, wei_op_t> >;
  // using incr_base_t = incrementalSolutionBase<
  //   romGalerkinImplicitJacobianPolicy, state_type, jacobian_type,
  //   model_type, time_type, sizer_type, basis_t>;
// private:
//   using incr_base_t::y0ptr_;
//   using incr_base_t::yFull_;

private:
  phi_op_t * phiOp_;
  wei_op_t * WOp_;
  jacobian_type JFOM_;
  state_type yFOM_;
  
public:
  romGalerkinImplicitJacobianPolicy(const state_type & y0fom,
				    const state_type & y0r,
				    phi_op_t & phiOp,
				    wei_op_t & WOp)
    : phiOp_(&phiOp), WOp_(&WOp),
      JFOM_(y0fom.size(), y0fom.size()), yFOM_(y0fom.size())
  {}

  ~romGalerkinImplicitJacobianPolicy() = default;
  
private:
  template <typename U = state_type,
	    typename T = jacobian_type,
	    typename std::enable_if<
	      core::meta::is_core_vector<U>::value==true &&
	      core::meta::is_coreMatrixWrapper<T>::value==true
	    >::type * = nullptr>
  void computeImpl(const U & y,
		   T & J,
		   model_type & model,
		   time_type t,
		   time_type dt)
  {
    // reconstruct the solution
    //    yFull_ = y;//*y0ptr_ + y;

    if (J.rows()==0)
      J.resize(y.size(), y.size());

    // reconstruct FOM state vector for computing space residual
    phiOp_->leftMultiply(y, yFOM_);

    // first eval space jacobian
    model.jacobian( *yFOM_.data(), *JFOM_.data(), t);

    WOp_->leftRightMultiply(JFOM_, J);

    // update from time discrete residual
    ode::impl::implicit_euler_time_discrete_jacobian(J, dt);
    
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
  }  

private:
  friend base_t;
  //  friend incr_base_t;

};//end class

  
}//end namespace exp
}//end namespace rom
#endif 
