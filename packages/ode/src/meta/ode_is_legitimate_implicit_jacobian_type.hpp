
#ifndef ODE_IS_LEGITIMATE_IMPLICIT_JACOBIAN_TYPE_HPP_
#define ODE_IS_LEGITIMATE_IMPLICIT_JACOBIAN_TYPE_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../containers/src/vector/containers_vector_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<typename jacobian_type, typename enable = void>
struct is_legitimate_jacobian_type : std::false_type{};

template<typename jacobian_type>
struct is_legitimate_jacobian_type<jacobian_type,
       typename std::enable_if<
	 containers::meta::is_matrix_wrapper<jacobian_type>::value or
	 containers::meta::is_multi_vector_wrapper<jacobian_type>::value
#ifdef HAVE_PYBIND11
	 or containers::meta::is_cstyle_array_pybind11<jacobian_type>::value
#endif
	 >::type
       > : std::true_type{};

}}} // namespace rompp::ode::meta
#endif
