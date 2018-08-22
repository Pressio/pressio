
#ifndef ROM_SVD_TRAITS_HPP_
#define ROM_SVD_TRAITS_HPP_

#include "rom_svd_forward_declarations.hpp"

namespace rom{
namespace details{

template<typename mat_type,
	 typename precond_type>
struct traits<<state_type,
					     ode_residual_type,
					     scalar_type,
					     model_type,
					     residual_policy_type>
	      >
{
  using state_t =  state_type;
  using ode_residual_t = ode_residual_type;
  using scalar_t = scalar_type;
  using model_t = model_type;
  using residual_policy_t = residual_policy_type;

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;

  using order_t = unsigned int;
  static constexpr order_t order_value = 1;    
};


////////////////////////////////////////////////////////////////

  
}//end namespace details
}//end namespace rom
#endif
