
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_NATIVE_EPETRA_VECTOR_META_HPP_
#define CONTAINERS_NATIVE_EPETRA_VECTOR_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

namespace rompp{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_epetra : std::false_type {};

template <typename T>
struct is_vector_epetra<T,
      typename
      std::enable_if<
	std::is_same<T,Epetra_Vector>::value
	>::type
      > : std::true_type{};

}}}//end namespace rompp::containers::meta
#endif
#endif
