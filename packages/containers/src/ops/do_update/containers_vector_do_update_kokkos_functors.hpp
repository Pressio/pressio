
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_CONTAINER_OPS_VECTOR_DO_UPDATE_KOKKOS_FUNCTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_VECTOR_DO_UPDATE_KOKKOS_FUNCTOR_HPP_

namespace pressio{ namespace containers{ namespace ops{ namespace impl{

template <typename T, typename sc_t>
struct DoUpdateTwoTermsFunctor {
  T v_ = {};
  T v1_ = {};
  T v2_ = {};
  sc_t a_ = ::pressio::utils::constants::zero<sc_t>();
  sc_t b_ = ::pressio::utils::constants::zero<sc_t>();
  sc_t c_ = ::pressio::utils::constants::zero<sc_t>();

  DoUpdateTwoTermsFunctor(T v, T v1, T v2,
			  sc_t a, sc_t b, sc_t c)
    : v_{v}, v1_{v1}, v2_{v2},
      a_{a}, b_{b}, c_{c}{}

  DoUpdateTwoTermsFunctor(T v, T v1, T v2,
			  sc_t b, sc_t c)
    : v_{v}, v1_{v1}, v2_{v2},
      b_{b}, c_{c}{}

  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
    v_(i) = a_*v_(i) + b_*v1_(i) + c_*v2_(i);
  }
};

}}}}//end namespace pressio::containers::ops::impl
#endif
#endif
