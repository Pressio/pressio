
#ifndef QR_BASE_SOLVE_BASE_HPP_
#define QR_BASE_SOLVE_BASE_HPP_

#include "../qr_meta.hpp"

namespace rompp{ namespace qr{

template<typename derived_t>
class QRSolveBase
  : private utils::details::CrtpBase<
  QRSolveBase<derived_t> >{

  using this_t = QRSolveBase<derived_t>;

  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_t>::type;

  friend utils::details::CrtpBase<this_t>;

public:
  template <
    typename vec_t,
    ::rompp::mpl::enable_if_t<
      containers::meta::is_vector_wrapper<vec_t>::value
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
