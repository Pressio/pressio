
#ifndef CORE_VECTOR_VECTOR_TRAITS_EXP_HPP_
#define CORE_VECTOR_VECTOR_TRAITS_EXP_HPP_

// #include <Eigen/Core>
// #include <type_traits>
// #include "Epetra_Vector.h"
// #include "core_ConfigDefs.hpp"
// #include "core_forward_declarations.hpp"

#include "../core_forward_declarations.hpp"
#include "../meta/core_vector_meta.hpp"
#include "../meta/core_meta_detect_typedefs.hpp"
#include "../meta/core_meta_detect_operators.hpp"


namespace core {
namespace details {


template <typename T, typename Enable = void>
struct vector_traits{

  using wrapped_type = T;

  static constexpr WrappedPackageIdentifier wrapped_package_identifier
      = WrappedPackageIdentifier::Undefined;
  
  static constexpr WrappedVectorIdentifier wrapped_vector_identifier
      = WrappedVectorIdentifier::Undefined;
  
  static constexpr bool is_dynamic = false;
  static constexpr int rows = -1;
};
//----------------------------------------------------

  
template <typename T>
struct vector_traits<
  core::Vector<T>,
  typename std::enable_if<
    std::is_base_of<
      Epetra_Vector,
      T
    >::value, void
  >::type
  >{

  using wrapped_type = T;

  static constexpr WrappedPackageIdentifier wrapped_package_identifier
      = WrappedPackageIdentifier::Trilinos;
  
  static constexpr WrappedVectorIdentifier wrapped_vector_identifier
      = WrappedVectorIdentifier::Epetra;
  
  static constexpr bool is_dynamic = true;
  static constexpr int rows = -1;

}; 
//----------------------------------------------------
  

template <typename T>
struct vector_traits<
  core::Vector<T>,
  typename std::enable_if<
    std::is_base_of<
      Eigen::Matrix<
        typename T::Scalar,
        T::RowsAtCompileTime,
        T::ColsAtCompileTime,
        T::Options,
        T::MaxRowsAtCompileTime,
        T::MaxColsAtCompileTime
      >, T
    >::value, void
  >::type
  >{
  
  using  wrapped_type = T;

  static constexpr WrappedPackageIdentifier wrapped_package_identifier
      = WrappedPackageIdentifier::Eigen;

  static constexpr WrappedVectorIdentifier wrapped_vector_identifier
      = WrappedVectorIdentifier::Eigen;
  
  static constexpr bool is_dynamic =
    T::RowsAtCompileTime == Eigen::Dynamic;

  static constexpr int rows =
    is_dynamic ? -1 : T::RowsAtCompileTime;
};
//----------------------------------------------------




} // end namespace details
} // end namespace core

#endif
