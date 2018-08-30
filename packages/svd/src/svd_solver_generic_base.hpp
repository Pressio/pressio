
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

public:

  void compute(matrix_t & mat, int nlsv, int nrsv){
    this->underlying().computeImpl(mat, nlsv, nrsv);
  }

  const leftSvec_t & cRefLeftSingularVectors() const {
    return this->underlying().cRefLeftSingularVectorsImpl();
  };

  const rightSvec_t & cRefRightSingularVectors() const {
    return this->underlying().cRefRightSingularVectorsImpl();
  };
  
  // const native_matrix_t & singularValues(){
  //   return this->underlying().singularValuesImpl();
  // };
    
private:    
  SolverBase() = default;
  ~SolverBase() = default;

private:
  friend derived_type;
  friend core::details::CrtpBase<SolverBase<derived_type>>;
  
};//end class

} // end namespace 
#endif
