
#ifndef CONTAINERS_OPS_META_HPP_
#define CONTAINERS_OPS_META_HPP_

#include "../vector/containers_vector_meta.hpp"
#include "../matrix/containers_matrix_meta.hpp"
#include "../multi_vector/containers_multi_vector_meta.hpp"


namespace pressio{ namespace containers{ namespace meta {


template <typename T1, typename T2, typename enable = void>
struct wrapper_pair_have_same_scalar : std::false_type {};

template <typename T1, typename T2>
struct wrapper_pair_have_same_scalar<T1,T2,
  ::pressio::mpl::enable_if_t<
    std::is_same<typename
		 containers::details::traits<T1>::scalar_t,
		 typename
		 containers::details::traits<T2>::scalar_t
		 >::value
    >
  > : std::true_type{};
//--------------------------------------------  


template <typename T1, typename T2,
	  typename T3, typename enable = void>
struct wrapper_triplet_have_same_scalar : std::false_type {};

template <typename T1, typename T2, typename T3>
struct wrapper_triplet_have_same_scalar<T1,T2,T3,
  ::pressio::mpl::enable_if_t<
    wrapper_pair_have_same_scalar<T1,T2>::value &&
    wrapper_pair_have_same_scalar<T2,T3>::value
    >
  > : std::true_type{};
//--------------------------------------------  

 
}}} // namespace pressio::containers::meta
#endif
