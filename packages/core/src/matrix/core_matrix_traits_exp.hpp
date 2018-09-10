
#ifndef CORE_MATRIX_MATRIX_TRAITS_EXP_HPP_
#define CORE_MATRIX_MATRIX_TRAITS_EXP_HPP_

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <type_traits>
#include "Epetra_RowMatrix.h"
#include "core_ConfigDefs.hpp"
#include "core_forward_declarations.hpp"

namespace core {
namespace details {

  
template <typename T, typename Enabled = void>
struct matrix_traits {
  typedef void wrapped_t;

  static constexpr WrappedPackageIdentifier wrapped_package_identifier
      = WrappedPackageIdentifier::Undefined;
  
  // static constexpr WrappedContainerIdentifier wrapped_container_identifier
  //     = WrappedContainerIdentifier::Undefined;

  static constexpr bool is_sparse = false;
};
//----------------------------------------------------


template <typename T>
struct matrix_traits<core::Matrix<T>,
  typename std::enable_if<
    std::is_base_of<
      Epetra_RowMatrix, T
    >::value, void
  >::type
> {

  static constexpr WrappedPackageIdentifier wrapped_package_identifier
      = WrappedPackageIdentifier::Trilinos;
  
  // static constexpr WrappedContainerIdentifier wrapped_container_identifier
  //     = WrappedContainerIdentifier::TrilinosEpetra;
  
  static constexpr bool is_sparse = true;
};
//----------------------------------------------------


template <typename T>
struct matrix_traits<
  core::Matrix<T>,
  typename std::enable_if<
    std::is_base_of<
      Eigen::SparseMatrix<
        typename T::Scalar,
        T::Options,
        typename T::StorageIndex
      >, T
    >::value, void 
  >::type
  >
{

  typedef T wrapped_type;

  static constexpr WrappedPackageIdentifier wrapped_package_identifier
      = WrappedPackageIdentifier::Eigen;

  // static constexpr WrappedContainerIdentifier wrapped_container_identifier
  //     = WrappedContainerIdentifier::Eigen;
  
  static constexpr bool is_sparse = true;
};
//----------------------------------------------------

  
} // end namespace details
} // end namespace core

#endif
