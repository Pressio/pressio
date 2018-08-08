
#ifndef ODE_AUX_DATA_HPP_
#define ODE_AUX_DATA_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace impl{

template<typename model_type,
	 typename residual_policy_type>
class ExpOdeAuxData
{
public:
  ExpOdeAuxData(model_type & mod,
		residual_policy_type & rpolo)
    : model_(&mod), residual_obj_(&rpolo){}
  ~ExpOdeAuxData() = default;

protected:
  model_type * model_;
  residual_policy_type * residual_obj_;
  
};
//---------------------------------------------

template<typename model_type,
	 typename scalar_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class ImpOdeAuxData
{

public:
  ImpOdeAuxData(model_type & mod,
		residual_policy_type & rpolo,
		jacobian_policy_type & jpolo)
    : model_(&mod), residual_obj_(&rpolo), jacobian_obj_(&jpolo){}
  ~ImpOdeAuxData() = default;

protected:
  model_type * model_;
  residual_policy_type * residual_obj_;
  jacobian_policy_type * jacobian_obj_;
  scalar_type t_;
  scalar_type dt_;
  
};

  
}//end namespace impl
}//end namespace ode  
  
#endif

  
