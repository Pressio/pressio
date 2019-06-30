
#ifndef CONTAINERS_CONTAINER_OPS_EIGEN_MAT_PRODUCT_EIGEN_VECTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_EIGEN_MAT_PRODUCT_EIGEN_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "../../matrix/containers_matrix_meta.hpp"
#include "../containers_eigen_ops_helper_impl.hpp"

namespace rompp{ namespace containers{ namespace ops{

/*---------------------------------------------------------
c = A b
- A is sparse matrix from eigen
- b is vector from eigen
- c is an eigen vector storing the result
---------------------------------------------------------*/

template <typename A_t, typename b_t, typename c_t,
 ::rompp::mpl::enable_if_t<
   containers::meta::is_sparse_matrix_wrapper_eigen<A_t>::value &&
   containers::meta::is_vector_wrapper_eigen<b_t>::value &&
   containers::meta::is_vector_wrapper_eigen<c_t>::value &&
   containers::meta::wrapper_pair_have_same_scalar<A_t,b_t>::value
   > * = nullptr
 >
void product(const A_t & A, const b_t & b, c_t & c){

  assert(A.cols() == b.size());
  assert(c.size() == A.rows());
  (*c.data()) = (*A.data()) * (*b.data());
}


template <typename A_t, typename b_t,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_sparse_matrix_wrapper_eigen<A_t>::value &&
    containers::meta::is_vector_wrapper_eigen<b_t>::value &&
    containers::meta::wrapper_pair_have_same_scalar<A_t, b_t>::value
    > * = nullptr
  >
auto product(const A_t & A, const b_t & b)
-> containers::Vector<
    Eigen::Matrix<typename containers::details::traits<A_t>::scalar_t,-1,1>
  >
{
  using sc_t = typename containers::details::traits<A_t>::scalar_t;
  containers::Vector<Eigen::Matrix<sc_t,-1,1>> c(A.rows());
  product(A,b,c);
  return c;
}


/*---------------------------------------------------------
c = A b
- A is dense matrix from eigen
- b is vector from eigen
---------------------------------------------------------*/

template <
  typename A_t, typename b_t, typename c_t,
  bool transposeA = false,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_dense_matrix_wrapper_eigen<A_t>::value and
    containers::meta::is_vector_wrapper_eigen<b_t>::value and
    containers::meta::is_vector_wrapper_eigen<c_t>::value and
    containers::meta::wrapper_triplet_have_same_scalar<A_t, b_t, c_t>::value and
    transposeA == false
    > * = nullptr
  >
void product(const A_t & A, const b_t & b, c_t & c){

  assert(A.cols() == b.size());
  assert(c.size() == A.rows());
  (*c.data()) = (*A.data()) * (*b.data());
}


template <
  typename A_t, typename b_t, typename c_t,
  bool transposeA = false,
  ::rompp::mpl::enable_if_t<
    containers::meta::is_dense_matrix_wrapper_eigen<A_t>::value and
    containers::meta::is_vector_wrapper_eigen<b_t>::value and
    containers::meta::is_vector_wrapper_eigen<c_t>::value and
    containers::meta::wrapper_triplet_have_same_scalar<A_t, b_t, c_t>::value and
    transposeA == true
    > * = nullptr
  >
void product(const A_t & A, const b_t & b, c_t & c){

  assert(A.rows() == b.size());
  assert(c.size() == A.cols());
  (*c.data()) = (*A.data()).transpose() * (*b.data());
}


template <typename A_t, typename b_t,
    ::rompp::mpl::enable_if_t<
      containers::meta::is_dense_matrix_wrapper_eigen<A_t>::value and
      containers::meta::is_vector_wrapper_eigen<b_t>::value and
      containers::meta::wrapper_pair_have_same_scalar<A_t,b_t>::value
      > * = nullptr
   >
auto product(const A_t & A, const b_t & b)
  -> containers::Vector<
    Eigen::Matrix<typename containers::details::traits<A_t>::scalar_t,-1,1>
    >
{
  using sc_t = typename containers::details::traits<A_t>::scalar_t;
  containers::Vector<Eigen::Matrix<sc_t,-1,1>> c(A.rows());
  product(A,b,c);
  return c;
}


}}}//end namespace rompp::containers;:ops
#endif
