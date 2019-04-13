
#ifndef SOLVERS_GET_MATRIX_SIZE_HELPER_HPP
#define SOLVERS_GET_MATRIX_SIZE_HELPER_HPP

#include "../../solvers_ConfigDefs.hpp"
// #include "../../solvers_system_traits.hpp"
// #include "../../solvers_meta_static_checks.hpp"
// #include "../helper_policies/solvers_converged_criterior_policy.hpp"
// #include "../helper_policies/solvers_norm_helper_policy.hpp"
// #include "../helper_policies/solvers_line_search_policy.hpp"

namespace rompp{ namespace solvers{ namespace impl{

template <typename T, typename enable = void>
struct MatrixGetSizeHelper;

template <typename T>
struct MatrixGetSizeHelper<
  T,
  ::rompp::mpl::enable_if_t<
    core::meta::is_core_multi_vector_wrapper<T>::value and
    core::details::traits<T>::is_shared_mem == false
    >
  >{
  static auto globalRows(const T & A) -> decltype(A.globalLength()){
    return A.globalLength();
  }
  static auto globalCols(const T & A) -> decltype(A.globalNumVectors()){
    return A.globalNumVectors();
  }
};


template <typename T>
struct MatrixGetSizeHelper<
  T,
  ::rompp::mpl::enable_if_t<
    core::meta::is_core_matrix_wrapper<T>::value and
    core::details::traits<T>::is_shared_mem == false
    >
  >{
  static auto globalRows(const T & A) -> decltype(A.globalRows()){
    return A.globalRows();
  }
  static auto globalCols(const T & A) -> decltype(A.globalCols()){
    return A.globalCols();
  }
};

template <typename T>
struct MatrixGetSizeHelper<
  T,
  ::rompp::mpl::enable_if_t<
    core::meta::is_core_matrix_wrapper<T>::value and
    core::details::traits<T>::is_shared_mem == true
    >
  >{
  static auto globalRows(const T & A) -> decltype(A.rows()){
    return A.rows();
  }
  static auto globalCols(const T & A) -> decltype(A.cols()){
    return A.cols();
  }
};

template <typename T>
struct MatrixGetSizeHelper<
  T,
  ::rompp::mpl::enable_if_t<
    core::meta::is_core_multi_vector_wrapper<T>::value and
    core::details::traits<T>::is_shared_mem == true
    >
  >{
  static auto globalRows(const T & A) -> decltype(A.length()){
    return A.length();
  }
  static auto globalCols(const T & A) -> decltype(A.numVectors()){
    return A.numVectors();
  }
};

}}} //end namespace rompp::solvers::impl
#endif
