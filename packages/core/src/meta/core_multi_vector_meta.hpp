
#ifndef CORE_MULTIVECTOR_META_MULTIVECTOR_META_HPP_
#define CORE_MULTIVECTOR_META_MULTIVECTOR_META_HPP_

#include "core_meta_basic.hpp"
#include "core_vector_meta.hpp"
#include "Epetra_MultiVector.h"

namespace core{
namespace meta {
 
template <typename T, typename enable = void>
struct is_multi_vector_epetra : std::false_type {};

template <typename T>
struct is_multi_vector_epetra<T,
  typename
   std::enable_if<
	std::is_same<T,Epetra_MultiVector>::value
   >::type
  > : std::true_type{};

//////////////////////
} // namespace meta
/////////////////////

#define STATIC_ASSERT_IS_MULTIVECTOR_EPETRA(TYPE) \
  static_assert( core::meta::is_multi_vector_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MULTIVECTOR_EPETRA")
#define STATIC_ASSERT_IS_NOT_MULTIVECTOR_EPETRA(TYPE) \
  static_assert( !core::meta::is_multi_vector_epetra<TYPE>::value, \
		 "THIS_IS_A_MULTIVECTOR_EPETRA")
  
  
/////////////////
} // namespace core
/////////////////

#endif
