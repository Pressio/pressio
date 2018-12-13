
#ifndef QR_SOLVER_BASE_HPP_
#define QR_SOLVER_BASE_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"

namespace rompp{ namespace qr{

template<typename derived_type,
	 typename matrix_type>
class QRSolverBase
  : private core::details::CrtpBase<
  QRSolverBase<derived_type, matrix_type>>{

  using this_t = QRSolverBase<derived_type, matrix_type>;
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

public:
  void compute(matrix_type & A){
    this->underlying().computeImpl(A);
  };

private:
  QRSolverBase() = default;
  ~QRSolverBase() = default;

};//end class

}}//end namespace rompp::qr
#endif
