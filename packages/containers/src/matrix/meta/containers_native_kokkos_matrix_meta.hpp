
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_NATIVE_KOKKOS_MATRIX_META_HPP_
#define CONTAINERS_NATIVE_KOKKOS_MATRIX_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include <KokkosSparse_CrsMatrix.hpp>

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_sparse_matrix_kokkos : std::false_type {};

template <typename T>
struct is_sparse_matrix_kokkos<
  T,
  typename
  std::enable_if<
    mpl::is_same<
      T,
      KokkosSparse::CrsMatrix<
	typename T::value_type,
	typename T::ordinal_type,
	typename T::execution_space,
	typename T::memory_traits,
	typename T::size_type
	>
      >::value
    >::type
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_dense_matrix_kokkos : std::false_type {};

template <typename T>
struct is_dense_matrix_kokkos<
  T,
  ::pressio::mpl::enable_if_t<
    // kokkos dense matrix is a view and has rank=2
    Kokkos::is_view<T>::value &&
    T::traits::rank==2
    >
  > : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif
