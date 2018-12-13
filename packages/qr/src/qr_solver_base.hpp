
#ifndef QR_SOLVER_BASE_HPP_
#define QR_SOLVER_BASE_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"

namespace rompp{ namespace qr{

template<typename derived_type,
	 typename Q_type,
	 typename R_type,
	 typename matrix_type>
class QRSolverBase
  : private core::details::CrtpBase<
  QRSolverBase<derived_type, Q_type,
	       R_type, matrix_type>>{

  using this_t = QRSolverBase<derived_type, Q_type,
			      R_type, matrix_type>;
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

public:
  void compute(matrix_type & A){
    this->underlying().computeImpl(A);
  };

  const Q_type & cRefQFactor() const {
    return this->underlying().cRefQFactorImpl();
    //return *Qmat_;
  }

  const R_type & cRefRFactor() const {
    return this->underlying().cRefRFactorImpl();
    //return *Rmat_;
  }

private:
  QRSolverBase() = default;
  ~QRSolverBase() = default;

};//end class

}}//end namespace rompp::qr
#endif
