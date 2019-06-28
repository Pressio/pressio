
#ifndef ALGEBRA_CONTAINER_OPS_MINMAX_MAX_VECTOR_HPP_
#define ALGEBRA_CONTAINER_OPS_MINMAX_MAX_VECTOR_HPP_

#include "../algebra_ops_meta.hpp"
#include "../../vector/algebra_vector_meta.hpp"

namespace rompp{ namespace algebra{ namespace ops{


//--------------------------------------------------------
//  eigen vector wrapper
//--------------------------------------------------------
template <typename vec_type,
  ::rompp::mpl::enable_if_t<
    ::rompp::algebra::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
auto max(const vec_type & a) -> typename details::traits<vec_type>::scalar_t
{
  return a.data()->maxCoeff();
}

// //--------------------------------------------------------
// //  blaze vector wrapper
// //--------------------------------------------------------
// #ifdef HAVE_BLAZE
// template <typename vec_type,
//   ::rompp::mpl::enable_if_t<
//     algebra::meta::is_dynamic_vector_wrapper_blaze<vec_type>::value or
//     algebra::meta::is_static_vector_wrapper_blaze<vec_type>::value
//     > * = nullptr
//   >
// auto max(const vec_type & a) -> typename details::traits<vec_type>::scalar_t
// {
// 	return blaze::max( *a.data() );
// }
// #endif// HAVE_BLAZE

// //--------------------------------------------------------
// //  armadillo vector wrapper
// //--------------------------------------------------------
// #ifdef HAVE_ARMADILLO
// template <typename vec_type,
//   ::rompp::mpl::enable_if_t<
//     algebra::meta::is_column_vector_wrapper_armadillo<vec_type>::value or
//     algebra::meta::is_row_vector_wrapper_armadillo<vec_type>::value
//     > * = nullptr
//   >
// auto max(const vec_type & a) -> typename details::traits<vec_type>::scalar_t
// {
// 	return arma::max( *a.data() );
// }
// #endif //HAVE_ARMADILLO

}}}//end namespace rompp::algebra::ops
#endif
