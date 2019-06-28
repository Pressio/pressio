
#ifndef ODE_IS_LEGITIMATE_IMPLICIT_JACOBIAN_TYPE_HPP_
#define ODE_IS_LEGITIMATE_IMPLICIT_JACOBIAN_TYPE_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../algebra/src/vector/algebra_vector_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<typename jacobian_type, typename enable = void>
struct is_legitimate_jacobian_type : std::false_type{};

template<typename jacobian_type>
struct is_legitimate_jacobian_type<jacobian_type,
       typename std::enable_if<
	 algebra::meta::is_algebra_matrix_wrapper<jacobian_type>::value or
	 algebra::meta::is_algebra_multi_vector_wrapper<jacobian_type>::value
#ifdef HAVE_PYBIND11
	 or algebra::meta::is_cstyle_array_pybind11<jacobian_type>::value
#endif
	 >::type
       > : std::true_type{};

}}} // namespace rompp::ode::meta
#endif
