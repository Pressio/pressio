
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_NATIVE_TPETRA_BLOCK_MULTI_VECTOR_META_HPP_
#define CONTAINERS_NATIVE_TPETRA_BLOCK_MULTI_VECTOR_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include <Tpetra_Experimental_BlockMultiVector_decl.hpp>

namespace rompp{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_tpetra_block : std::false_type {};

template <typename T>
struct is_multi_vector_tpetra_block<
  T,
  ::rompp::mpl::enable_if_t<
    std::is_same<T,
		 Tpetra::Experimental::BlockMultiVector<
		   typename T::impl_scalar_type,
		   typename T::local_ordinal_type,
		   typename T::global_ordinal_type,
		   typename T::node_type
		   >
		 >::value
    >
  > : std::true_type{};


}}}//end namespace rompp::containers::meta
#endif
#endif
