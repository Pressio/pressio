
#ifndef solvers_has_create_hessian_method_return_results_hpp_
#define solvers_has_create_hessian_method_return_results_hpp_

namespace pressio{ namespace solvers{ namespace meta {

template<typename T, typename h_t, typename enable = void>
struct has_const_create_hessian_method_return_result : std::false_type{};

template<typename T, typename h_t>
struct has_const_create_hessian_method_return_result
<T, h_t,
 ::pressio::mpl::enable_if_t<
   ::pressio::mpl::is_same<
     h_t,
     decltype( std::declval<T const>().createHessian() )
     >::value
   >
 > : std::true_type{};

}}} // namespace pressio::solvers::meta
#endif
