
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_NATIVE_TEUCHOS_VECTOR_META_HPP_
#define CONTAINERS_NATIVE_TEUCHOS_VECTOR_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include "Teuchos_SerialDenseVector.hpp"

namespace rompp{ namespace containers{ namespace meta {

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

}}}//end namespace rompp::containers::meta
#endif
#endif
