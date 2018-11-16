
#ifndef ODE_AUX_DATA_HPP_
#define ODE_AUX_DATA_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename model_type,
	 typename residual_policy_type>
class ExpOdeAuxData
{
public:
  ExpOdeAuxData(const model_type & mod,
		const residual_policy_type & rpolo)
    : model_(&mod), residual_obj_(&rpolo){}

  ExpOdeAuxData(const residual_policy_type & rpolo)
    : model_(nullptr), residual_obj_(&rpolo){}

  ~ExpOdeAuxData() = default;

protected:
  const model_type * model_ = nullptr;
  const residual_policy_type * residual_obj_ = nullptr;

};
//---------------------------------------------

template<typename model_type,
	 typename scalar_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class ImpOdeAuxData
{

public:
  ImpOdeAuxData(const model_type & mod,
		const residual_policy_type & rpolo,
		const jacobian_policy_type & jpolo)
    : model_(&mod), residual_obj_(&rpolo), jacobian_obj_(&jpolo){}
  ~ImpOdeAuxData() = default;

protected:
  const model_type * model_ = nullptr;
  const residual_policy_type * residual_obj_ = nullptr;
  const jacobian_policy_type * jacobian_obj_ = nullptr;
  scalar_type t_;
  scalar_type dt_;
};


}}}//end namespace rompp::ode::impl
#endif
