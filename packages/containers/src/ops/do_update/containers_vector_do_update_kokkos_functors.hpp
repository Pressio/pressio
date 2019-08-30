
#ifdef HAVE_KOKKOS
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

template <typename T, typename sc_t>
struct DoUpdateThreeTermsFunctor {
  T v_ = {};
  T v1_ = {};
  T v2_ = {};
  T v3_ = {};
  sc_t a_ = ::pressio::utils::constants::zero<sc_t>();
  sc_t b_ = ::pressio::utils::constants::zero<sc_t>();
  sc_t c_ = ::pressio::utils::constants::zero<sc_t>();
  sc_t d_ = ::pressio::utils::constants::zero<sc_t>();

  DoUpdateThreeTermsFunctor(T v, T v1, T v2, T v3,
			    sc_t a, sc_t b, sc_t c, sc_t d)
    : v_{v}, v1_{v1}, v2_{v2}, v3_{v3},
      a_{a}, b_{b}, c_{c}, d_{d}{}

  DoUpdateThreeTermsFunctor(T v, T v1, T v2, T v3,
			    sc_t b, sc_t c, sc_t d)
    : v_{v}, v1_{v1}, v2_{v2}, v3_{v3},
      b_{b}, c_{c}, d_{d}{}

  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
    v_(i) = a_*v_(i) + b_*v1_(i) + c_*v2_(i) + d_*v3_(i);
  }
};



template <typename T, typename sc_t>
struct DoUpdateFourTermsFunctor {
  T v_ = {};
  T v1_ = {};
  T v2_ = {};
  T v3_ = {};
  T v4_ = {};
  sc_t a_ = ::pressio::utils::constants::zero<sc_t>();
  sc_t b_ = ::pressio::utils::constants::zero<sc_t>();
  sc_t c_ = ::pressio::utils::constants::zero<sc_t>();
  sc_t d_ = ::pressio::utils::constants::zero<sc_t>();
  sc_t e_ = ::pressio::utils::constants::zero<sc_t>();

  DoUpdateFourTermsFunctor(T v, T v1, T v2, T v3, T v4,
			   sc_t a, sc_t b, sc_t c, sc_t d, sc_t e)
    : v_{v}, v1_{v1}, v2_{v2}, v3_{v3}, v4_{v4},
      a_{a}, b_{b}, c_{c}, d_{d}, e_{e}{}

  DoUpdateFourTermsFunctor(T v, T v1, T v2, T v3, T v4,
			   sc_t b, sc_t c, sc_t d, sc_t e)
    : v_{v}, v1_{v1}, v2_{v2}, v3_{v3}, v4_{v4},
      b_{b}, c_{c}, d_{d}, e_{e}{}

  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
    v_(i) = a_*v_(i) + b_*v1_(i) + c_*v2_(i) + d_*v3_(i) + e_*v4_(i);
  }
};


}}}}//end namespace pressio::containers::ops::impl
#endif
#endif
