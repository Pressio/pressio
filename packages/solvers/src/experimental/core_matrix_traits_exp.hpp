
#ifndef CORE_MATRIX_TRAITS_EXP_HPP_
#define CORE_MATRIX_TRAITS_EXP_HPP_

#include <Eigen/Core>


namespace core {
namespace details {

template <typename T, typename Enabled = void>
struct matrix_traits {
  constexpr bool is_eigen = false;
  constexpr bool is_trilinos = false;
};

template <T>
struct matrix_traits<core::Matrix<T>,
  std::enable_if<std::is_base_of<EPetra_RowMatrix, T>, void>::type> {
  constexpr bool is_eigen = false;
  constexpr bool is_trilinos = true;
};

template <T>
struct matrix_traits<core::Matrix<T>,
  std::enable_if<std::is_base_of<Eigen::SparseMatrix<typename T::Scalar,
      T::Options,
      typename T::StorageIndex
    >, T
  >, void>::type> {
  constexpr bool is_eigen = true;
  constexpr bool is_trilinos = false;
};
	
} // end namespace details
} // end namespace core

#endif
