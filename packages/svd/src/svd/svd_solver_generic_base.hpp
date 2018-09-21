
#ifndef SVD_SOLVER_GENERIC_BASE_HPP_
#define SVD_SOLVER_GENERIC_BASE_HPP_

#include "svd_solver_traits.hpp"

namespace svd{
    
template<typename derived_type>
class SolverBase
  : private core::details::CrtpBase<SolverBase<derived_type>>
{

private:
  using matrix_t = typename svd::details::svd_traits<derived_type>::matrix_t;
  using sc_t = typename svd::details::svd_traits<derived_type>::scalar_t;
  using leftSvec_t = typename svd::details::svd_traits<derived_type>::lsv_t;
  using rightSvec_t = typename svd::details::svd_traits<derived_type>::rsv_t;
  using sval_t = typename svd::details::svd_traits<derived_type>::sval_t;

public:

  template<svdType envalue,
	   typename std::enable_if<
	     envalue==svdType::truncated
	     >::type * = nullptr>
  void compute(matrix_t & mat, sc_t tol, int t){
    this->underlying().template computeImpl<envalue>(mat, tol, t);
  }

  const leftSvec_t & cRefLeftSingularVectors() const {
    return this->underlying().cRefLeftSingularVectorsImpl();
  };

  const rightSvec_t & cRefRightSingularVectors() const {
    return this->underlying().cRefRightSingularVectorsImpl();
  };
  
  const sval_t & singularValues() const{
    return this->underlying().singularValuesImpl();
  };
    
private:    
  SolverBase() = default;
  ~SolverBase() = default;

private:
  friend derived_type;
  friend core::details::CrtpBase<SolverBase<derived_type>>;
  
};//end class

} // end namespace 
#endif
