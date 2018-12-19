
#ifndef QR_SOLVER_BASE_HPP_
#define QR_SOLVER_BASE_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "qr_rfactor_solve_impl.hpp"

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
  void computeThin(matrix_type & A){
    this->underlying().computeThinImpl(A);
  }

  template <typename vector_in_t,
	    typename vector_out_t>
  void project(const vector_in_t & vecIn,
	       vector_out_t & vecOut) const{
    this->underlying().projectImpl(vecIn, vecOut);
  }

  template <typename vector_t>
  void solve(const vector_t & rhs, vector_t & y){
    this->underlying().solveImpl(rhs, y);
  }

  const Q_type & cRefQFactor() const {
    return this->underlying().cRefQFactorImpl();
  }

  const R_type & cRefRFactor() const {
    return this->underlying().cRefRFactorImpl();
  }

private:
  QRSolverBase() = default;
  ~QRSolverBase() = default;

  std::shared_ptr<Q_type> Qmat_ = nullptr;
  std::shared_ptr<R_type> Rmat_ = nullptr;

};//end class

}}//end namespace rompp::qr
#endif
