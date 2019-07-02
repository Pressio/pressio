
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_NATIVE_TPETRA_VECTOR_META_HPP_
#define CONTAINERS_NATIVE_TPETRA_VECTOR_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include "Tpetra_Vector.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_tpetra : std::false_type {};

template <typename T>
struct is_vector_tpetra<T,
      typename
      std::enable_if<
	std::is_same<T,
		     Tpetra::Vector<
		       typename T::impl_scalar_type,
		       typename T::local_ordinal_type,
		       typename T::global_ordinal_type,
		       typename T::node_type
		       >
		     >::value
	>::type
      > : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif
