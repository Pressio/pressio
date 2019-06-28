
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_NATIVE_TEUCHOS_VECTOR_META_HPP_
#define ALGEBRA_NATIVE_TEUCHOS_VECTOR_META_HPP_

#include "../../meta/algebra_meta_basic.hpp"
#include "Teuchos_SerialDenseVector.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_vector_teuchos : std::false_type {};

template <typename T>
struct is_dense_vector_teuchos<T,
      ::rompp::mpl::enable_if_t<
	std::is_same<T,
	  Teuchos::SerialDenseVector<typename T::ordinalType,
				     typename T::scalarType>
	  >::value
	>
      > : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
#endif
