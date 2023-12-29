#ifndef EXPRESSIONS_TEST_HELPERS_HPP_
#define EXPRESSIONS_TEST_HELPERS_HPP_

#include "pressio/expressions.hpp"

namespace helpers
{

template <typename T>
void fill_matrix(T & A)
{
  A(0,0) = 1.2; A(0,1) = 2.;  A(0,2) = 3.;   A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.2; A(1,2) = 7.;   A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.2; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.;  A(3,3) = 16.;
}


template<class T>
void fill_matrix_seq(T & A){
  A(0,0) = 1.;  A(0,1) = 2.;  A(0,2) = 3.;  A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.;  A(1,2) = 7.;  A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.; A(3,3) = 16.;
  A(4,0) = 17.; A(4,1) = 18.; A(4,2) = 19.; A(4,3) = 20.;
  A(5,0) = 21.; A(5,1) = 22.; A(5,2) = 23.; A(5,3) = 24.;
}

template <typename Scalar, typename T>
void check_span_traits(T obj)
{
  using expr_t = decltype(pressio::span(obj, 0, 1));

  static_assert(pressio::Traits<expr_t>::rank == 1);
  static_assert(std::is_same_v<typename pressio::Traits<expr_t>::scalar_type, Scalar>);
  static_assert(std::is_same_v<typename pressio::Traits<expr_t>::reference_type, Scalar &>);
}

template <typename Scalar, typename T>
void check_subspan_traits(T obj)
{
  using pair_t = std::pair<std::size_t, std::size_t>;
  using expr_t = decltype(pressio::subspan(obj, pair_t{0, 1}, pair_t{0,1}));

  static_assert(pressio::Traits<expr_t>::rank == 2);
  static_assert(std::is_same_v<typename pressio::Traits<expr_t>::scalar_type, Scalar>);
  static_assert(std::is_same_v<typename pressio::Traits<expr_t>::reference_type, Scalar &>);
}

} // namespace helpers

#endif  // EXPRESSIONS_TEST_HELPERS_HPP_
