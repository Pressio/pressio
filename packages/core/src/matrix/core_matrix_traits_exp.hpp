
#ifndef CORE_MATRIX_TRAITS_EXP_HPP_
#define CORE_MATRIX_TRAITS_EXP_HPP_

#include <Eigen/Core>
#include <type_traits>

#include "Epetra_RowMatrix.h"

#include "core_forward_declarations.hpp"


namespace core {
namespace details {

template <typename T, typename Enabled = void>
struct matrix_traits {
  typedef void wrapped_t;
  static constexpr bool is_eigen = false;
  static constexpr bool is_trilinos = false;
  static constexpr bool is_sparse = false;
};


template <typename T>
struct matrix_traits<core::matrix<T>,
  typename std::enable_if<
    std::is_base_of<
      Epetra_RowMatrix, T
    >::value, void
  >::type
> {
  static constexpr bool is_eigen = false;
  static constexpr bool is_trilinos = true;
  static constexpr bool is_sparse = true;
};


template <typename T>
struct matrix_traits<
  core::matrix<T>,
  typename std::enable_if<
    std::is_base_of<
      Eigen::SparseMatrix<
        typename T::Scalar,
        T::Options,
        typename T::StorageIndex
      >, T
    >::value, void 
  >::type
> {
  typedef T wrapped_t;
  static constexpr bool is_eigen = true;
  static constexpr bool is_trilinos = false;
  static constexpr bool is_sparse = true;
};
	
} // end namespace details
} // end namespace core

#endif