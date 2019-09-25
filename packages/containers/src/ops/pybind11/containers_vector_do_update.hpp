
#ifdef HAVE_PYBIND11
#ifndef CONTAINERS_SRC_OPS_PYBIND11_VECTOR_DO_UPDATE_HPP_
#define CONTAINERS_SRC_OPS_PYBIND11_VECTOR_DO_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace pressio{ namespace containers{ namespace ops{


//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, scalar_t a,
	       const T & v1, scalar_t b){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = a*v.at(i) + b*v1.at(i);
  }
}


template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, const T & v1, const scalar_t b){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = b*v1.at(i);
  }
}



//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){
  std::cout << " do_update pybind11 two terms " << std::endl;
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }
  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = a*v.at(i) + b*v1.at(i) + c*v2.at(i);
  }
}


template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = b*v1.at(i) + c*v2.at(i);
  }
}



//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size() and vsz != v3.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = a*v.at(i) + b*v1.at(i) + c*v2.at(i) + d*v3.at(i);
  }
}


template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size() and vsz != v3.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = b*v1.at(i) + c*v2.at(i) + d*v3.at(i);
  }
}



//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size()
      and vsz != v3.size() and vsz != v4.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = a*v.at(i) + b*v1.at(i) + c*v2.at(i) + d*v3.at(i) + e*v4.at(i);
  }
}


template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size()
      and vsz != v3.size() and vsz != v4.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = b*v1.at(i) + c*v2.at(i) + d*v3.at(i) + e*v4.at(i);
  }
}


}}}//end namespace pressio::containers::ops
#endif
#endif
