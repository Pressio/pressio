
#ifndef CORE_EIGEN_MAT_MAT_IMPL_HPP_
#define CORE_EIGEN_MAT_MAT_IMPL_HPP_

#include "../core_ops_meta.hpp"
#include "../../matrix/core_matrix_meta.hpp"

namespace rompp{ namespace core{ namespace ops{ namespace impl{


template <typename TA, typename TB, typename TC, typename enable = void>
struct are_matrices_admissible_for_eigen_mat_mat : std::false_type{};

template <typename TA, typename TB, typename TC>
struct are_matrices_admissible_for_eigen_mat_mat<
  TA, TB, TC,
  core::meta::enable_if_t<
    core::meta::is_eigen_matrix_wrapper<TA>::value and
    core::meta::is_eigen_matrix_wrapper<TB>::value and
    core::meta::is_eigen_matrix_wrapper<TC>::value and
    core::meta::wrapper_triplet_have_same_scalar<TA,TB,TC>::value and
    // for now just enable for dynamic
    core::details::traits<TA>::is_dynamic and
    core::details::traits<TB>::is_dynamic and
    core::details::traits<TC>::is_dynamic
    >
  > : std::true_type{};
//----------------------------------------------------------------

template <typename TA, typename TB>
struct are_matrices_admissible_for_eigen_mat_mat<
  TA, TB, void,
  core::meta::enable_if_t<
    core::meta::is_eigen_matrix_wrapper<TA>::value and
    core::meta::is_eigen_matrix_wrapper<TB>::value and
    core::meta::wrapper_pair_have_same_scalar<TA,TB>::value and
    // for now just enable for dynamic
    core::details::traits<TA>::is_dynamic and
    core::details::traits<TB>::is_dynamic
    >
  > : std::true_type{};
//----------------------------------------------------------------


template <typename TA, typename TB, typename enable = void>
struct returnTypeHelper;

// TA = dense, TB = dense, their product is dense
template <typename TA>
struct returnTypeHelper<
  TA, TA,
  core::meta::enable_if_t<
    core::meta::is_eigen_dense_matrix_wrapper<TA>::value
    >
  >{
  using prod_type = TA;
};

// TA = sparse, TB = sparse, their product is dense
template <typename TA>
struct returnTypeHelper<
  TA, TA,
  core::meta::enable_if_t<
    core::meta::is_eigen_sparse_matrix_wrapper<TA>::value
    >
  >{
  using prod_type = TA;
};

// TA = sparse, TB = dense, their product is dense
template <typename TA, typename TB>
struct returnTypeHelper<
  TA, TB,
  core::meta::enable_if_t<
    core::meta::is_eigen_sparse_matrix_wrapper<TA>::value and
    core::meta::is_eigen_dense_matrix_wrapper<TB>::value
    >
  >{
  using prod_type = TB;
};

// TA = dense, TB = sparse, their product is dense
template <typename TA, typename TB>
struct returnTypeHelper<
  TA, TB,
  core::meta::enable_if_t<
    core::meta::is_eigen_dense_matrix_wrapper<TA>::value and
    core::meta::is_eigen_sparse_matrix_wrapper<TB>::value
    >
  >{
  using prod_type = TA;
};
//-------------------------------------------------------------


template <bool transpose_A, bool transpose_B>
struct eig_mat_mat_product;

// Default is C = A * B
template <>
struct eig_mat_mat_product<false, false>{

  template <typename TA, typename TB, typename TC>
  void operator()(const TA & A, const TB & B, TC & C) const{
    assert(A.cols() == B.rows());
    assert(C.rows() == A.rows());
    assert(C.cols() == B.cols());
    (*C.data()) = (*A.data()) * (*B.data());
  }

  template <typename TA, typename TB>
  auto operator()(const TA & A, const TB & B) const
    -> typename returnTypeHelper<TA, TB>::prod_type
  {
    using prod_ret_t = typename returnTypeHelper<TA, TB>::prod_type;
    prod_ret_t C(A.rows(),B.cols());
    (*this)(A,B,C);
    return C;
  }

};


// C = A^T * B
template <>
struct eig_mat_mat_product<true, false>{

  template <typename TA, typename TB, typename TC >
  void operator()(const TA & A, const TB & B, TC & C) const{
    assert(A.rows() == B.rows());
    assert(C.rows() == A.cols());
    assert(C.cols() == B.cols());
    (*C.data()) = (*A.data()).transpose() * (*B.data());
  }

  template <typename TA, typename TB>
  auto operator()(const TA & A, const TB & B) const
    -> typename returnTypeHelper<TA, TB>::prod_type
  {
    using prod_ret_t = typename returnTypeHelper<TA, TB>::prod_type;
    prod_ret_t C(A.cols(),B.cols());
    (*this)(A,B,C);
    return C;
  }
};


// C = A * B^T
template <>
struct eig_mat_mat_product<false, true>{

  template <typename TA, typename TB, typename TC >
  void operator()(const TA & A, const TB & B, TC & C) const{
    assert(A.cols() == B.cols());
    assert(C.rows() == A.rows());
    assert(C.cols() == B.rows());
    (*C.data()) = (*A.data()) * (*B.data()).transpose();
  }

  template <typename TA, typename TB>
  auto operator()(const TA & A, const TB & B) const
    -> typename returnTypeHelper<TA, TB>::prod_type
  {
    using prod_ret_t = typename returnTypeHelper<TA, TB>::prod_type;
    prod_ret_t C(A.rows(), B.rows());
    (*this)(A,B,C);
    return C;
  }
};


}}}}//end namespace rompp::core::ops::impl
#endif
