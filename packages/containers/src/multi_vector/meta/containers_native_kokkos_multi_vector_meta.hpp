
#ifdef HAVE_KOKKOS
#ifndef CONTAINERS_NATIVE_KOKKOS_MULTI_VECTOR_META_HPP_
#define CONTAINERS_NATIVE_KOKKOS_MULTI_VECTOR_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
// #include <Tpetra_MultiVector_decl.hpp>

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_kokkos : std::false_type {};

template <typename T>
struct is_multi_vector_kokkos<
  T,
  ::pressio::mpl::enable_if_t<
    // kokkos MV is a view and has rank=2
    Kokkos::is_view<T>::value &&
    T::traits::rank==2
    >
  > : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif
