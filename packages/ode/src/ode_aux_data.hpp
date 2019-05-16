
#ifndef ODE_AUX_DATA_HPP_
#define ODE_AUX_DATA_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename model_type,
	 typename residual_policy_type>
class ExpOdeAuxData
{
protected:
  ExpOdeAuxData(const model_type & mod,
		const residual_policy_type & rpolo)
    : model_(&mod), residual_obj_(&rpolo){}

  ExpOdeAuxData(const residual_policy_type & rpolo)
    : model_(nullptr), residual_obj_(&rpolo){}

  ~ExpOdeAuxData() = default;

  const model_type * model_		     = nullptr;
  const residual_policy_type * residual_obj_ = nullptr;

};
//---------------------------------------------

template<
  typename model_type,
  typename scalar_type
  >
struct ImpOdeAuxData{

  const model_type & model_;
  scalar_type t_	= {};
  scalar_type dt_	= {};

  ImpOdeAuxData(const model_type & model)
    : model_(model){}

  ~ImpOdeAuxData() = default;

};


}}}//end namespace rompp::ode::impl
#endif
