
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_NATIVE_EPETRA_MULTI_VECTOR_META_HPP_
#define ALGEBRA_NATIVE_EPETRA_MULTI_VECTOR_META_HPP_

#include "../../meta/algebra_meta_basic.hpp"
#include "Epetra_MultiVector.h"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_epetra : std::false_type {};

template <typename T>
struct is_multi_vector_epetra<T,
  typename
   std::enable_if<
    std::is_same<T,Epetra_MultiVector>::value
   >::type
  > : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
#endif
