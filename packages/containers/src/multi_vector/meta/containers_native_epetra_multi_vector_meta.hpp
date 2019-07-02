
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_NATIVE_EPETRA_MULTI_VECTOR_META_HPP_
#define CONTAINERS_NATIVE_EPETRA_MULTI_VECTOR_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include "Epetra_MultiVector.h"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_epetra : std::false_type {};

template <typename T>
struct is_multi_vector_epetra<T,
  typename
   std::enable_if<
    std::is_same<T,Epetra_MultiVector>::value
   >::type
  > : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif
