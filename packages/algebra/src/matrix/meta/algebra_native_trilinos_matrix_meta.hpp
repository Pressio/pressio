
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_NATIVE_TRILINOS_MATRIX_META_HPP_
#define ALGEBRA_NATIVE_TRILINOS_MATRIX_META_HPP_

#include "../../meta/algebra_meta_basic.hpp"
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_LocalMap.h>
#include <Tpetra_CrsMatrix_decl.hpp>
#include "Teuchos_SerialDenseMatrix.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_sparse_matrix_tpetra : std::false_type {};

template <typename T>
struct is_sparse_matrix_tpetra<T,
      typename
      std::enable_if<
  std::is_same<T,
         Tpetra::CrsMatrix<
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
struct is_sparse_matrix_epetra
  : std::false_type {};

template<typename T>
struct is_sparse_matrix_epetra<T,
    typename std::enable_if<
      std::is_same<T, Epetra_CrsMatrix>::value
      >::type >
  : std::true_type{};

//-------------------------------------------------

template <typename T, typename enable = void>
struct is_dense_matrix_epetra
  : std::false_type {};

template<typename T>
struct is_dense_matrix_epetra<T,
    typename std::enable_if<
      std::is_same<T, Epetra_MultiVector>::value
      >::type >
  : std::true_type{};
//-------------------------------------------------

template <typename T, typename enable = void>
struct is_dense_matrix_teuchos : std::false_type {};

template <typename T>
struct is_dense_matrix_teuchos<T,
    typename std::enable_if<
	std::is_same<T,
	  Teuchos::SerialDenseMatrix<typename T::ordinalType,
				     typename T::scalarType>
	  >::value
	>::type
      > : std::true_type{};
//-------------------------------------------------

template <typename T, typename enable = void>
struct is_dense_matrix_teuchos_rcp : std::false_type {};

template <typename T>
struct is_dense_matrix_teuchos_rcp<T,
    typename std::enable_if<
      is_teuchos_rcp<T>::value and
      is_dense_matrix_teuchos<typename T::element_type>::value
	>::type
      > : std::true_type{};


}}}//end namespace rompp::algebra::meta
#endif
#endif
