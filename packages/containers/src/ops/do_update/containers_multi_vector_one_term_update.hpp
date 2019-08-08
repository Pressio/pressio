
#ifndef CONTAINERS_CONTAINER_OPS_MULTI_VECTOR_ONE_TERM_UPDATE_HPP_
#define CONTAINERS_CONTAINER_OPS_MULTI_VECTOR_ONE_TERM_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#ifdef HAVE_KOKKOS
#include<KokkosBlas1_axpby.hpp>
#endif

//----------------------------------------------------------------------
//  overloads for computing:
//	MV = a * MV + b * MV1
// where MV is a multivector wrapper
//----------------------------------------------------------------------

namespace pressio{ namespace containers{ namespace ops{

//----------------------------------------------------------------------
// enable for Eigen wrappers
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
  for (auto i=0; i<mv.length(); i++){
    for (auto j=0; j<mv.numVectors(); j++)
      mv(i,j) = a*mv(i,j) + b*mv1(i,j);
  }
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper<T>::value
    > * = nullptr
  >
void do_update(T & mv, const T & mv1, const scalar_t & b)
{
  assert( mv.numVectors() == mv1.numVectors() );
  assert( mv.lenght() == mv1.length() );
  for (auto i=0; i<mv.length(); i++){
    for (auto j=0; j<mv.numVectors(); j++)
      mv(i,j) = b*mv1(i,j);
  }
}

//---------------------------------------------------------------------
// enable for kokkos wrapper
//---------------------------------------------------------------------
#ifdef HAVE_KOKKOS
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & mv, const scalar_t &a,
	       const T & mv1, const scalar_t &b)
{
  KokkosBlas::axpby(b, *mv1.data(), a, *mv.data());
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & mv, const T & mv1, const scalar_t &b)
{
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
  KokkosBlas::axpby(b, *mv1.data(), zero, *mv.data());
}
#endif


}}}//end namespace pressio::containers::ops
#endif
