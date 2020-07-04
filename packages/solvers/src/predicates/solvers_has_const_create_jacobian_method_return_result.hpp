
#ifndef solvers_has_create_jacobian_method_return_results_hpp_
#define solvers_has_create_jacobian_method_return_results_hpp_

namespace pressio{ namespace solvers{ namespace meta {

template<typename T, typename jac_t, typename enable = void>
struct has_const_create_jacobian_method_return_result : std::false_type{};

template<typename T, typename jac_t>
struct has_const_create_jacobian_method_return_result
<T, jac_t,
 ::pressio::mpl::enable_if_t<
   ::pressio::mpl::is_same<
     jac_t,
     decltype( std::declval<T const>().createJacobian() )
     >::value
   >
 > : std::true_type{};

}}} // namespace pressio::solvers::meta
#endif
