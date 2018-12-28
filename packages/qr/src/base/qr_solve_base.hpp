
#ifndef QR_BASE_SOLVE_BASE_HPP_
#define QR_BASE_SOLVE_BASE_HPP_

#include "../qr_meta.hpp"

namespace rompp{ namespace qr{

template<typename derived_t>
class QRSolveBase
  : private core::details::CrtpBase<
  QRSolveBase<derived_t> >{

  using this_t = QRSolveBase<derived_t>;
  friend derived_t;
  friend core::details::CrtpBase<this_t>;

public:
  template <
    typename vec_t,
    core::meta::enable_if_t<
      core::meta::is_core_vector_wrapper<vec_t>::value
      >* = nullptr
    >
  void solve(const vec_t & rhs, vec_t & y)const {
    this->underlying().solveImpl(rhs, y);
  }

private:
  QRSolveBase() = default;
  ~QRSolveBase() = default;

};

}}//end namespace rompp::qr
#endif
