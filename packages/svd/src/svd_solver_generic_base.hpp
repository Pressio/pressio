
#ifndef SVD_SOLVER_GENERIC_BASE_HPP_
#define SVD_SOLVER_GENERIC_BASE_HPP_

#include "svd_solver_traits.hpp"

namespace svd{
    
template<typename derived_type>
class SolverBase
  : private core::details::CrtpBase<SolverBase<derived_type>>
{

private:
  using sc_t = typename svd::details::traits<derived_type>::scalar_t;
  using wrap_solver_t = typename svd::details::traits<derived_type>::wrapped_solver_t;
  //using native_matrix_t = typename svd::details::traits<derived_type>::native_matrix_t;
  // using u_matrix_type = typename svd::details::traits<derived_type>::u_matrix_t;
  //using v_matrix_type = typename svd::details::traits<derived_type>::v_matrix_t;

public:

  template <typename matrix_in_type>
  void compute(const matrix_in_type & mat,
	       bool isHermitian, 
	       int nev){
    this->underlying().computeImpl(mat);
  };

  // const native_matrix_t & singularValues(){
  //   return this->underlying().singularValuesImpl();
  // };

  // const u_matrix_type & leftSingularVectors() const {
  //   return this->underlying().leftSingularVectorsImpl();
  // };

  // const v_matrix_type & rightSingularVectors() const {
  //   return this->underlying().rightSingularVectorsImpl();
  // };
  
  // wrap_solver_t const & getConstRefToWrappedSolver() const {
  //   return this->underlying().getConstRefToWrappedSolverImpl();
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
