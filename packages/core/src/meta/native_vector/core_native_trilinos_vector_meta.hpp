
#ifdef HAVE_TRILINOS
#ifndef CORE_NATIVE_TRILINOS_VECTOR_META_HPP_
#define CORE_NATIVE_TRILINOS_VECTOR_META_HPP_

#include "../core_meta_basic.hpp"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Tpetra_Vector.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_SerialDenseVector.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_epetra : std::false_type {};

template <typename T>
struct is_vector_epetra<T,
      typename
      std::enable_if<
	std::is_same<T,Epetra_Vector>::value
	>::type
      > : std::true_type{};
//--------------------------------------------

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
//--------------------------------------------

template <typename T, typename enable = void>
struct is_dense_vector_teuchos : std::false_type {};

template <typename T>
struct is_dense_vector_teuchos<T,
      core::meta::enable_if_t<
	std::is_same<T,
	  Teuchos::SerialDenseVector<typename T::ordinalType,
				     typename T::scalarType>
	  >::value
	>
      > : std::true_type{};
//--------------------------------------------

template <typename T, typename enable = void>
struct is_vector_kokkos : std::false_type {};

template <typename T>
struct is_vector_kokkos<T,
	 core::meta::enable_if_t<
	   // kokkos vector is it is a view and has rank=1
	   Kokkos::is_view<T>::value &&
	   T::traits::rank==1>
      > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
#endif
