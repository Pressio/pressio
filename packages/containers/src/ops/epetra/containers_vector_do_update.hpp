
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_SRC_OPS_EPETRA_VECTOR_DO_UPDATE_HPP_
#define CONTAINERS_SRC_OPS_EPETRA_VECTOR_DO_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t a, const T & v1, const scalar_t b){
  v = a*v + b*v1;
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const T & v1, const scalar_t  b){
  v = b*v1;
}


//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){
  v = a*v + b*v1 + c*v2;
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){
  v = b*v1 + c*v2;
}


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d){
  v = a*v + b*v1 + c*v2 + d*v3;
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d){
  v = b*v1 + c*v2 + d*v3;
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){
  v = a*v + b*v1 + c*v2 + d*v3 + e*v4;
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){
  v = b*v1 + c*v2 + d*v3 + e*v4;
}


}}}//end namespace pressio::containers::ops
#endif
#endif
