#ifndef CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DO_UPDATE_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DO_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

//----------------------------------------------------------------------
//  overloads for computing: MV = a * MV + b * MV1
// where MV is an eigen multivector wrapper
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<T>::value
    > * = nullptr
  >
void do_update(T & mv, const scalar_t &a,
	       const T & mv1, const scalar_t &b)
{
  assert( mv.numVectors() == mv1.numVectors() );
  assert( mv.length() == mv1.length() );
  for (decltype(mv.length()) i=0; i<mv.length(); i++){
    for (decltype(mv.numVectors()) j=0; j<mv.numVectors(); j++)
      mv(i,j) = a*mv(i,j) + b*mv1(i,j);
  }
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<T>::value
    > * = nullptr
  >
void do_update(T & mv, const T & mv1, const scalar_t & b)
{
  assert( mv.numVectors() == mv1.numVectors() );
  assert( mv.lenght() == mv1.length() );
  for (decltype(mv.length()) i=0; i<mv.length(); i++){
    for (decltype(mv.numVectors()) j=0; j<mv.numVectors(); j++)
      mv(i,j) = b*mv1(i,j);
  }
}

}}}//end namespace pressio::containers::ops
#endif
