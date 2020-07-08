
#ifndef solvers_has_create_gradient_method_return_results_hpp_
#define solvers_has_create_gradient_method_return_results_hpp_

namespace pressio{ namespace solvers{ namespace predicates {

template<typename T, typename g_t, typename enable = void>
struct has_const_create_gradient_method_return_result : std::false_type{};

template<typename T, typename g_t>
struct has_const_create_gradient_method_return_result
<T, g_t,
 ::pressio::mpl::enable_if_t<
   ::pressio::mpl::is_same<
     g_t,
     decltype( std::declval<T const>().createGradient() )
     >::value
   >
 > : std::true_type{};

}}} // namespace pressio::solvers::predicates
#endif
