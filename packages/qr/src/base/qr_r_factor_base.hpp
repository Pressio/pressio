
#ifndef QR_BASE_R_FACTOR_BASE_HPP_
#define QR_BASE_R_FACTOR_BASE_HPP_

#include "../qr_meta.hpp"

namespace rompp{ namespace qr{

template<typename derived_t,
	 typename R_t>
class RFactorBase
  : private core::details::CrtpBase<
  RFactorBase<derived_t, R_t>>{

  using this_t = RFactorBase<derived_t, R_t>;
  friend derived_t;
  friend core::details::CrtpBase<this_t>;

public:
  const R_t & cRefRFactor() const {
    return this->underlying().cRefRFactorImpl();
  }

private:
  RFactorBase() = default;
  ~RFactorBase() = default;

  mutable std::shared_ptr<R_t> Rmat_ = nullptr;

};


}}//end namespace rompp::qr
#endif
