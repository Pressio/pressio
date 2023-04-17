#ifndef OPS_FIXTURES_SHARED_LEVEL2_HPP_
#define OPS_FIXTURES_SHARED_LEVEL2_HPP_

// Note: all structures must be accessible from host
template <
  typename TransMode,
  typename AType,
  typename XType,
  typename YType,
  typename ScalarType
  >
void vanilla_gemv(TransMode /* trans */, ScalarType alpha,
                  AType A, XType x, ScalarType beta, YType &y) {
  const bool trans = std::is_same<TransMode, ::pressio::transpose>::value;
  const auto x_size = ::pressio::ops::extent(x, 0);
  const auto y_size = ::pressio::ops::extent(y, 0);
  assert(y_size == ::pressio::ops::extent(A, trans ? 1 : 0));
  assert(x_size == ::pressio::ops::extent(A, trans ? 0 : 1));
  if (beta == 0.) ::pressio::ops::set_zero(y);
  else ::pressio::ops::scale(y, beta);
  if (alpha == 0.) return;
  for (size_t i = 0; i < y_size; ++i) {
    for (size_t j = 0; j < x_size; ++j) {
      const auto a = trans ? A(j, i) : A(i, j);
      y(i) += alpha * a * x(j);
    }
  }
}

#endif // OPS_FIXTURES_SHARED_LEVEL2_HPP_
