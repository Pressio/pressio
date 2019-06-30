
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_NATIVE_KOKKOS_VECTOR_META_HPP_
#define CONTAINERS_NATIVE_KOKKOS_VECTOR_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include "Kokkos_Core.hpp"

namespace rompp{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_kokkos : std::false_type {};

template <typename T>
struct is_vector_kokkos<T,
	 ::rompp::mpl::enable_if_t<
	   // kokkos vector is it is a view and has rank=1
	   Kokkos::is_view<T>::value &&
	   T::traits::rank==1>
      > : std::true_type{};

}}}//end namespace rompp::containers::meta
#endif
#endif
