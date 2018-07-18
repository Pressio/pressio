
#ifndef CORE_VECTOR_TRAITS_EXP_HPP_
#define CORE_VECTOR_TRAITS_EXP_HPP_

#include <Eigen/Core>
#include <type_traits>

#include "Epetra_Vector.h"

#include "core_ConfigDefs.hpp"
#include "core_forward_declarations.hpp"


namespace core {
namespace details {


template <typename T, typename Enable = void>
struct vector_traits {
  typedef T wrapped_type;
  static constexpr WrappedClass vector_class = WrappedClass::Undefined;
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
  static constexpr WrappedClass vector_class = WrappedClass::Trilinos;
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
  static constexpr WrappedClass vector_class = WrappedClass::Eigen;
  static constexpr bool is_dynamic = T::RowsAtCompileTime == Eigen::Dynamic;
  static constexpr int rows = is_dynamic ? -1 : T::RowsAtCompileTime;
};


template <typename T, typename U>
struct same_vector_structure {
  static constexpr bool value = vector_traits<T>::is_dynamic || vector_traits<U>::is_dynamic || vector_traits<T>::rows == vector_traits<U>::rows;
};


template <typename T, typename U>
struct are_vector_compatible {
  static constexpr bool valid_vector = vector_traits<T>::vector_class != WrappedClass::Undefined;
  static constexpr bool same_type = vector_traits<T>::vector_class == vector_traits<U>::vector_class;
  static constexpr bool same_structure = vector_traits<T>::is_dynamic || vector_traits<U>::is_dynamic || vector_traits<T>::rows == vector_traits<U>::rows;
  static constexpr bool value = valid_vector && same_type && same_structure;
};


} // end namespace details
} // end namespace core

#endif