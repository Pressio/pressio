
#ifndef ODE_IS_LEGITIMATE_IMPLICIT_JACOBIAN_TYPE_HPP_
#define ODE_IS_LEGITIMATE_IMPLICIT_JACOBIAN_TYPE_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../core/src/vector/core_vector_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<typename jacobian_type, typename enable = void>
struct is_legitimate_jacobian_type : std::false_type{};

template<typename jacobian_type>
struct is_legitimate_jacobian_type<jacobian_type,
       typename std::enable_if<
	 core::meta::is_core_matrix_wrapper<jacobian_type>::value or
	 core::meta::is_core_multi_vector_wrapper<jacobian_type>::value
	 >::type
       > : std::true_type{};

}}} // namespace rompp::ode::meta
#endif
