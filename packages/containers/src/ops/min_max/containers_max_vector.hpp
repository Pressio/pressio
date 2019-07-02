
#ifndef CONTAINERS_CONTAINER_OPS_MINMAX_MAX_VECTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_MINMAX_MAX_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{


//--------------------------------------------------------
//  eigen vector wrapper
//--------------------------------------------------------
template <typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_eigen<vec_type>::value
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
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_dynamic_vector_wrapper_blaze<vec_type>::value or
//     containers::meta::is_static_vector_wrapper_blaze<vec_type>::value
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
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_column_vector_wrapper_armadillo<vec_type>::value or
//     containers::meta::is_row_vector_wrapper_armadillo<vec_type>::value
//     > * = nullptr
//   >
// auto max(const vec_type & a) -> typename details::traits<vec_type>::scalar_t
// {
// 	return arma::max( *a.data() );
// }
// #endif //HAVE_ARMADILLO

}}}//end namespace pressio::containers::ops
#endif
