
#ifndef ODE_EULER_STANDARD_POLICIES_HPP_
#define ODE_EULER_STANDARD_POLICIES_HPP_

#include "ode_ConfigDefs.hpp"

#include "../impl/ode_euler_implicit_impl.hpp"


namespace ode{  
namespace policies{  

  
template <typename model_type>
class policyBase{
public:
  policyBase(model_type * model) : model_(model){}  
protected:
   model_type * model_;
};

  
template <typename derived_type, typename model_type>
class explicitEulerResidualPolicyBase : protected policyBase<model_type>
{
public:
  explicitEulerResidualPolicyBase(model_type * model)
    : policyBase<model_type>(model){}

  template<typename state_type, typename residual_type>
  void compute(const state_type & y, residual_type & R){
    this->policy()->computeImpl(y,R);
  } 
private:
  derived_type * policy(){
    return static_cast< derived_type* >( this );
  }
  const derived_type * policy() const{
    return static_cast< const derived_type* >( this );
  }
};


  
template <typename derived_type, typename model_type>
class implicitEulerResidualPolicyBase : private policyBase<model_type>
{
public:
  implicitEulerResidualPolicyBase(model_type * model)
    : policyBase<model_type>(model){}

  template<typename state_type,
	   typename residual_type,
	   typename time_type>
  void compute(const state_type & y, const state_type & ynm1,
	       residual_type & R, time_type dt){
    this->policy()->computeImpl(y,ynm1,R,dt);
  } 
private:
  derived_type * policy(){
    return static_cast< derived_type* >( this );
  }
  const derived_type * policy() const{
    return static_cast< const derived_type* >( this );
  }
};


template <typename derived_type, typename model_type>
class implicitEulerJacobianPolicyBase : private policyBase<model_type>
{
public:
  implicitEulerJacobianPolicyBase(model_type * model)
    : policyBase<model_type>(model){}

  template<typename state_type, typename jacobian_type, typename time_type>
  void compute(const state_type & y, jacobian_type & J, time_type dt){
    this->policy()->computeImpl(y,J,dt);
  } 
private:
  derived_type * policy(){
    return static_cast< derived_type* >( this );
  }
  const derived_type * policy() const{
    return static_cast< const derived_type* >( this );
  }
};
///===================================================



template <typename model_type>
class explicitEulerStandardResidual
  : public explicitEulerResidualPolicyBase< explicitEulerStandardResidual<model_type>,
					    model_type>
{
public:
  explicitEulerStandardResidual(model_type * model)
    : explicitEulerResidualPolicyBase<explicitEulerStandardResidual<model_type>,
				      model_type>(model){}

  template<typename state_type, typename residual_type>
  void computeImpl(const state_type & y, residual_type & R){
    this->model_->residual(y,R);
  }
};


template <typename model_type>
class implicitEulerStandardResidual
  : public implicitEulerResidualPolicyBase<implicitEulerStandardResidual<model_type>,
					   model_type>
{
public:
  implicitEulerStandardResidual(model_type * model)
    : implicitEulerResidualPolicyBase<implicitEulerStandardResidual<model_type>,
				      model_type>(model){}

  template<typename state_type, typename residual_type,typename time_type>
  void compute(const state_type & y, const state_type & ynm1, residual_type & R, time_type dt)
  {
    // first eval RHS
    this->model_->residual(y,R);
    // then fix residual based on time stepping features
    ode::impl::implicit_euler_residual_impl(y, ynm1, R, dt);
  }
};


template <typename model_type>
class implicitEulerStandardJacobian
  : public implicitEulerJacobianPolicyBase<implicitEulerStandardJacobian<model_type>,
					   model_type>
{
public:
  implicitEulerStandardJacobian(model_type * model)
    : implicitEulerJacobianPolicyBase<implicitEulerStandardJacobian<model_type>,
				      model_type>(model){}

  template<typename state_type, typename jacobian_type, typename time_type>
  void compute(const state_type & y, jacobian_type & J, time_type dt)
  {
    // first eval jac
    this->model_->jacobian(y,J);
    // then fix it based on time stepping features
    ode::impl::implicit_euler_jacobian_impl(y, J, dt);
  }
};
    



}//end namespace polices
}//end namespace ode  
#endif 

