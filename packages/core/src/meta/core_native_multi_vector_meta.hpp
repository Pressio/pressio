
#ifndef CORE_NATIVE_MULTIVECTOR_META_MULTIVECTOR_META_HPP_
#define CORE_NATIVE_MULTIVECTOR_META_MULTIVECTOR_META_HPP_

#include "core_meta_basic.hpp"
#include "core_native_vector_meta.hpp"
// #include "Epetra_MultiVector.h"

namespace rompp{
namespace core{
namespace meta {
 
#ifdef HAVE_TRILINOS
template <typename T, typename enable = void>
struct is_multi_vector_epetra : std::false_type {};

template <typename T>
struct is_multi_vector_epetra<T,
  typename
   std::enable_if<
	std::is_same<T,Epetra_MultiVector>::value
   >::type
  > : std::true_type{};
#endif


} // namespace meta
} // namespace core

}//end namespace rompp
#endif
