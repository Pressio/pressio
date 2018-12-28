
#ifndef QR_BASE_OUT_OF_PLACE_BASE_HPP_
#define QR_BASE_OUT_OF_PLACE_BASE_HPP_

#include "../qr_meta.hpp"

namespace rompp{ namespace qr{


template<typename derived_t,
	 typename matrix_t,
	 typename Q_t>
class QROutOfPlaceBase
  : private core::details::CrtpBase<
  QROutOfPlaceBase<derived_t, matrix_t, Q_t>>{

  using this_t = QROutOfPlaceBase<derived_t, matrix_t, Q_t>;
  friend derived_t;
  friend core::details::CrtpBase<this_t>;

public:
  void computeThin(matrix_t & A){
    this->underlying().computeThinImpl(A);
  }

  const Q_t & cRefQFactor() const {
    return this->underlying().cRefQFactorImpl();
  }

  template <
    typename vec_in_t,
    typename vec_out_t,
    core::meta::enable_if_t<
      core::meta::is_core_vector_wrapper<vec_in_t>::value and
      core::meta::is_core_vector_wrapper<vec_out_t>::value and
      meta::is_legitimate_vector_type_for_qr_project<vec_in_t,
						     Q_t>::value
      >* = nullptr
    >
  void project(const vec_in_t & vecIn,
	       vec_out_t & vecOut) const{
    this->underlying().projectImpl(vecIn, vecOut);
  }

private:
  QROutOfPlaceBase() = default;
  ~QROutOfPlaceBase() = default;

};

}}//end namespace rompp::qr
#endif
