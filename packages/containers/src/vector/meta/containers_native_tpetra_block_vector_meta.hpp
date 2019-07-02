
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_NATIVE_TPETRA_BLOCK_VECTOR_META_HPP_
#define CONTAINERS_NATIVE_TPETRA_BLOCK_VECTOR_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include <Tpetra_Experimental_BlockVector_decl.hpp>

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_tpetra_block : std::false_type {};

template <typename T>
struct is_vector_tpetra_block<T,
      typename
      std::enable_if<
	std::is_same<T,
		     Tpetra::Experimental::BlockVector<
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
