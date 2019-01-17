
#ifndef CORE_CONTAINER_OPS_MINMAX_MAX_VECTOR_HPP_
#define CORE_CONTAINER_OPS_MINMAX_MAX_VECTOR_HPP_

#include "../core_ops_meta.hpp"
#include "../../vector/core_vector_meta.hpp"

namespace rompp{ namespace core{ namespace ops{


//--------------------------------------------------------
//  eigen vector wrapper
//--------------------------------------------------------
template <typename vec_type,
  core::meta::enable_if_t<
    core::meta::is_eigen_vector_wrapper<vec_type>::value
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
//   core::meta::enable_if_t<
//     core::meta::is_blaze_dynamic_vector_wrapper<vec_type>::value or
//     core::meta::is_blaze_static_vector_wrapper<vec_type>::value
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
//   core::meta::enable_if_t<
//     core::meta::is_column_vector_armadillo_wrapper<vec_type>::value or
//     core::meta::is_row_vector_armadillo_wrapper<vec_type>::value
//     > * = nullptr
//   >
// auto max(const vec_type & a) -> typename details::traits<vec_type>::scalar_t
// {
// 	return arma::max( *a.data() );
// }
// #endif //HAVE_ARMADILLO

}}}//end namespace rompp::core::ops
#endif
