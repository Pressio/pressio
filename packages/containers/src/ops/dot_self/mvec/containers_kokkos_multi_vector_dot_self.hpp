
#ifdef HAVE_KOKKOS
#ifndef CONTAINERS_KOKKOS_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_KOKKOS_MULTI_VECTOR_DOT_SELF_HPP_

#include "../../containers_ops_meta.hpp"
#include "../../../multi_vector/containers_multi_vector_meta.hpp"
#include "KokkosBlas3_gemm.hpp"

namespace pressio{ namespace containers{ namespace ops{

// kokkos multivector wrapper (A) dot self
// this is equivalent to doing A^T * A
// returns a dense container kokkos wrapper matrix
// in the form of a containers::Matrix< Kokkos::View<> >

template <
  typename mvec_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & A,
	      containers::Matrix<
		Kokkos::View<
		typename containers::details::traits<mvec_t>::scalar_t**,
		typename containers::details::traits<mvec_t>::layout,
		typename containers::details::traits<mvec_t>::execution_space
	       >
	      > & C)
{

  using sc_t = typename containers::details::traits<mvec_t>::scalar_t;
  constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
  constexpr auto one = ::pressio::utils::constants::one<sc_t>();
  const char ctA = 'T';
  const char ctB = 'N';
  KokkosBlas::gemm(&ctA, &ctB, one, *A.data(), *A.data(), zero, *C.data());
}


template <typename mvec_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_t>::value
    > * = nullptr
  >
auto dot_self(const mvec_t & mvA)
  -> containers::Matrix<
    Kokkos::View<
      typename containers::details::traits<mvec_t>::scalar_t **,
      typename containers::details::traits<mvec_t>::layout,
      typename containers::details::traits<mvec_t>::execution_space
      >
    >{

  using sc_t	  = typename containers::details::traits<mvec_t>::scalar_t;
  using layout    = typename containers::details::traits<mvec_t>::layout;
  using exe_space = typename containers::details::traits<mvec_t>::execution_space;

  using dm_t = Kokkos::View<sc_t**, layout, exe_space>;
  using res_t = containers::Matrix<dm_t>;

  const auto numVecsA = mvA.numVectors();
  res_t C("dot_self_mat_res", numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}

}}} // end namespace pressio::containers::ops
#endif
#endif
