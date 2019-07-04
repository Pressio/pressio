
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_KOKKOS_MULTI_VECTOR_PRODUCT_KOKKOS_VECTOR_HPP_
#define CONTAINERS_KOKKOS_MULTI_VECTOR_PRODUCT_KOKKOS_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace containers{ namespace ops{

template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_type>::value and
    containers::meta::is_vector_wrapper_kokkos<vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & A,
	 const vec_type & b,
	 containers::Vector<
	   Kokkos::View<
	     typename containers::details::traits<mvec_type>::scalar_t*,
	     typename containers::details::traits<mvec_type>::layout,
	     typename containers::details::traits<mvec_type>::execution_space
	   >
	 > & c)
{
  static_assert(meta::kokkos_wrapper_pair_have_same_exe_space<mvec_type, vec_type>::value,
		"product: MV and vec types need to have same execution space" );

  assert( A.data()->extent(1) == b.data()->extent(0) );
  assert( c.data()->extent(0) == A.data()->extent(0) );

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
  constexpr auto one = ::pressio::utils::constants::one<sc_t>();

  const char ctA = 'N';
  KokkosBlas::gemv(&ctA, one, *A.data(), *b.data(), zero, *c.data());
}


// result is built and returned
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_type>::value and
    containers::meta::is_vector_wrapper_kokkos<vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
auto product(const mvec_type & mvA, const vec_type & vecB)
-> containers::Vector<
    Kokkos::View<
      typename containers::details::traits<mvec_type>::scalar_t*,
      typename containers::details::traits<mvec_type>::layout,
      typename containers::details::traits<mvec_type>::execution_space
      >
  >{

  static_assert(meta::kokkos_wrapper_pair_have_same_exe_space<mvec_type, vec_type>::value,
		"product: MV and vec types need to have same execution space" );

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  using layout    = typename containers::details::traits<mvec_type>::layout;
  using exe_space = typename containers::details::traits<mvec_type>::execution_space;

  using v_t = Kokkos::View<sc_t*, layout, exe_space>;
  using res_t = containers::Vector<v_t>;

  res_t c("product_res", mvA.data()->extent(0));
  assert( mvA.data()->extent(1) == vecB.data()->extent(0) );
  product(mvA, vecB, c);
  return c;
}
//--------------------------------------------------------

}}} // end namespace pressio::containers::ops
#endif
#endif
