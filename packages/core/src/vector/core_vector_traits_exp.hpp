
#ifndef CORE_VECTOR_TRAITS_EXP_HPP_
#define CORE_VECTOR_TRAITS_EXP_HPP_

#include <Eigen/Core>
#include <type_traits>

#include "Epetra_Vector.h"

#include "core_forward_declarations.hpp"


namespace core {
namespace details {

template <typename T, typename Enable = void>
struct vector_traits {
  typedef T wrapped_type;
  static constexpr bool is_eigen = false;
  static constexpr bool is_trilinos = false;
  static constexpr bool is_dynamic = false;
  static constexpr int rows = -1;
};


template <typename T>
struct vector_traits<
  core::vector<T>,
  typename std::enable_if<
    std::is_base_of<
      Epetra_Vector,
      T
    >::value, void
  >::type
> {
  typedef T wrapped_type;
  static constexpr bool is_eigen = false;
  static constexpr bool is_trilinos = true;
  static constexpr bool is_dynamic = true;
  static constexpr int rows = -1;
}; 
  

template <typename T>
struct vector_traits<
  core::vector<T>,
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
> {
  typedef T wrapped_type;
  static constexpr bool is_eigen = true;
  static constexpr bool is_trilinos = false;
  static constexpr bool is_dynamic = T::RowsAtCompileTime == Eigen::Dynamic;
  static constexpr int rows = is_dynamic ? -1 : T::RowsAtCompileTime;
};

} // end namespace details
} // end namespace core

#endif