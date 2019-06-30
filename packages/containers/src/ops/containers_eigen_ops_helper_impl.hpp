
#ifndef CONTAINERS_EIGEN_OPS_HELPER_IMPL_HPP_
#define CONTAINERS_EIGEN_OPS_HELPER_IMPL_HPP_

#include "containers_ops_meta.hpp"
#include "../matrix/containers_matrix_meta.hpp"

namespace rompp{ namespace containers{ namespace ops{ namespace impl{


template <typename TA, typename TB, typename enable = void>
struct eigenMatMatProdRetTypeHelper;

// TA = dense, TB = dense, their product is dense
template <typename TA>
struct eigenMatMatProdRetTypeHelper<
  TA, TA,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_dense_matrix_wrapper_eigen<TA>::value
    >
  >{
  using prod_type = TA;
};

// TA = sparse, TB = sparse, their product is sparse
template <typename TA>
struct eigenMatMatProdRetTypeHelper<
  TA, TA,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_sparse_matrix_wrapper_eigen<TA>::value
    >
  >{
  using prod_type = TA;
};

// TA = sparse, TB = dense, their product is dense
template <typename TA, typename TB>
struct eigenMatMatProdRetTypeHelper<
  TA, TB,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_sparse_matrix_wrapper_eigen<TA>::value and
    containers::meta::is_dense_matrix_wrapper_eigen<TB>::value
    >
  >{
  using prod_type = TB;
};

// TA = dense, TB = sparse, their product is dense
template <typename TA, typename TB>
struct eigenMatMatProdRetTypeHelper<
  TA, TB,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_dense_matrix_wrapper_eigen<TA>::value and
    containers::meta::is_sparse_matrix_wrapper_eigen<TB>::value
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
    assert(A.data()->cols() == B.data()->rows());
    assert(C.data()->rows() == A.data()->rows());
    assert(C.data()->cols() == B.data()->cols());
    (*C.data()) = (*A.data()) * (*B.data());
  }

  template <typename TA, typename TB>
  auto operator()(const TA & A, const TB & B) const
    -> typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type
  {
    using prod_ret_t = typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type;
    prod_ret_t C(A.data()->rows(),B.data()->cols());
    (*this)(A,B,C);
    return C;
  }

};


// C = A^T * B
template <>
struct eig_mat_mat_product<true, false>{

  template <typename TA, typename TB, typename TC >
  void operator()(const TA & A, const TB & B, TC & C) const{
    assert(A.data()->rows() == B.data()->rows());
    assert(C.data()->rows() == A.data()->cols());
    assert(C.data()->cols() == B.data()->cols());
    (*C.data()) = (*A.data()).transpose() * (*B.data());
  }

  template <typename TA, typename TB>
  auto operator()(const TA & A, const TB & B) const
    -> typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type
  {
    using prod_ret_t = typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type;
    prod_ret_t C(A.data()->cols(),B.data()->cols());
    (*this)(A,B,C);
    return C;
  }
};


// C = A * B^T
template <>
struct eig_mat_mat_product<false, true>{

  template <typename TA, typename TB, typename TC >
  void operator()(const TA & A, const TB & B, TC & C) const{
    assert(A.data()->cols() == B.data()->cols());
    assert(C.data()->rows() == A.data()->rows());
    assert(C.data()->cols() == B.data()->rows());
    (*C.data()) = (*A.data()) * (*B.data()).transpose();
  }

  template <typename TA, typename TB>
  auto operator()(const TA & A, const TB & B) const
    -> typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type
  {
    using prod_ret_t = typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type;
    prod_ret_t C(A.data()->rows(), B.data()->rows());
    (*this)(A,B,C);
    return C;
  }
};


}}}}//end namespace rompp::containers::ops::impl
#endif
