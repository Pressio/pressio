
#ifndef CORE_OPS_META_HPP_
#define CORE_OPS_META_HPP_

#include "../matrix/core_matrix_meta.hpp"
#include "../vector/core_vector_meta.hpp"
#include "../multi_vector/core_multi_vector_meta.hpp"


namespace rompp{
namespace core{
namespace meta {


template <typename T1, typename T2, typename enable = void>
struct wrapper_pair_have_same_scalar : std::false_type {};

template <typename T1, typename T2>
struct wrapper_pair_have_same_scalar<T1,T2,
  core::meta::enable_if_t<
    std::is_same<typename
		 core::details::traits<T1>::scalar_t,
		 typename
		 core::details::traits<T2>::scalar_t
		 >::value
    >
  > : std::true_type{};
//--------------------------------------------  


template <typename T1, typename T2,
	  typename T3, typename enable = void>
struct wrapper_triplet_have_same_scalar : std::false_type {};

template <typename T1, typename T2, typename T3>
struct wrapper_triplet_have_same_scalar<T1,T2,T3,
  core::meta::enable_if_t<
    wrapper_pair_have_same_scalar<T1,T2>::value &&
    wrapper_pair_have_same_scalar<T2,T3>::value
    >
  > : std::true_type{};
//--------------------------------------------  

  
 
} // namespace meta
} // namespace core

}//end namespace rompp
#endif
