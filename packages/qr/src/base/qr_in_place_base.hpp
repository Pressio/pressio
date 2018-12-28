
#ifndef QR_BASE_IN_PLACE_BASE_HPP_
#define QR_BASE_IN_PLACE_BASE_HPP_

#include "../qr_meta.hpp"

namespace rompp{ namespace qr{

template<typename derived_t,
	 typename matrix_t>
class QRInPlaceBase
  : private core::details::CrtpBase<
  QRInPlaceBase<derived_t, matrix_t>>{

  using this_t = QRInPlaceBase<derived_t, matrix_t>;
  friend derived_t;
  friend core::details::CrtpBase<this_t>;

public:
  void computeThin(matrix_t & A){
    this->underlying().computeThinImpl(A);
  }

private:
  QRInPlaceBase() = default;
  ~QRInPlaceBase() = default;

};

}}//end namespace rompp::qr
#endif
