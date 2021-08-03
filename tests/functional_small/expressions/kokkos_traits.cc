
#include <gtest/gtest.h>
#include "pressio_expressions.hpp"

TEST(expressions_kokkos, span_traits)
{
  {
    using T = Kokkos::View<double*>;
    T a("a", 10);
    using expr_t = decltype(pressio::span(a, 0, 1));
    static_assert(pressio::Traits<expr_t>::is_static,"");

    using cT = Kokkos::View<const double*>;
    using ref_t = typename cT::reference_type;
    static_assert(std::is_const<typename std::remove_reference<ref_t>::type >::value,"");

    cT b = a; 
    using expr1_t = decltype(pressio::span(b, 0, 1));
    static_assert(pressio::Traits<expr1_t>::is_static,"");
    using ref2_t = typename pressio::Traits<expr1_t>::reference_type;
    static_assert(std::is_const<typename std::remove_reference<ref2_t>::type>::value,"");
  }

  {
    using T = Kokkos::View<double*>;
    T a0("a", 10);
    const T a = a0;
    using expr_t = decltype(pressio::span(a, 0, 1));
    static_assert(pressio::Traits<expr_t>::is_static,"");
  }

  {
    using T = Kokkos::View<double[4]>;
    T a("a");
    using expr_t = decltype(pressio::span(a, 0, 1));
    static_assert(pressio::Traits<expr_t>::is_static,"");
  }

  {
    using T = Kokkos::View<double[4]>;
    T a0("a");
    const T a = a0;
    using expr_t = decltype(pressio::span(a, 0, 1));
    static_assert(pressio::Traits<expr_t>::is_static,"");
  }
}

TEST(expressions_kokkos, subspan_traits)
{
  {
    using T = Kokkos::View<double**>;
    T a("a", 10, 10);
    using pair_t = std::pair<int,int>;
    using expr_t = decltype(pressio::subspan(a, pair_t{0, 1}, pair_t{0,1}));
    static_assert(pressio::Traits<expr_t>::is_static,"");

    using cT = Kokkos::View<const double**>;
    using ref_t = typename cT::reference_type;
    static_assert(std::is_const<typename std::remove_reference<ref_t>::type >::value,"");

    cT b = a; 
    using expr1_t = decltype(pressio::subspan(b, pair_t{0, 1}, pair_t{0,1}));
    static_assert(pressio::Traits<expr1_t>::is_static,"");
    using ref2_t = typename pressio::Traits<expr1_t>::reference_type;
    static_assert(std::is_const<typename std::remove_reference<ref2_t>::type>::value,"");
  }

  {
    using T = Kokkos::View<double**>;
    const T a("a", 10, 10);
    using pair_t = std::pair<int,int>;
    using expr_t = decltype(pressio::subspan(a, pair_t{0, 1}, pair_t{0,1}));
    static_assert(pressio::Traits<expr_t>::is_static,"");
  }

  {
    using T = Kokkos::View<double*[4]>;
    T a("a", 10);
    using pair_t = std::pair<int,int>;
    using expr_t = decltype(pressio::subspan(a, pair_t{0, 1}, pair_t{0,1}));
    static_assert(pressio::Traits<expr_t>::is_static,"");
  }
}


TEST(expressions_kokkos, diag_traits)
{
  {
    using T = Kokkos::View<double**>;
    T a("a", 10, 10);
    using expr_t = decltype(pressio::diag(a));
    static_assert(pressio::Traits<expr_t>::is_static,"");

    using cT = Kokkos::View<const double**>;
    using ref_t = typename cT::reference_type;
    static_assert(std::is_const<typename std::remove_reference<ref_t>::type >::value,"");

    cT b = a; 
    using expr1_t = decltype(pressio::diag(b));
    static_assert(pressio::Traits<expr1_t>::is_static,"");
    using ref2_t = typename pressio::Traits<expr1_t>::reference_type;
    static_assert(std::is_const<typename std::remove_reference<ref2_t>::type>::value,"");
  }

  {
    using T = Kokkos::View<double**>;
    const T a("a", 10, 10);
    using expr_t = decltype(pressio::diag(a));
    static_assert(pressio::Traits<expr_t>::is_static,"");
  }

  {
    using T = Kokkos::View<double*[10]>;
    T a("a", 10);
    using expr_t = decltype(pressio::diag(a));
    static_assert(pressio::Traits<expr_t>::is_static,"");
  }
}

