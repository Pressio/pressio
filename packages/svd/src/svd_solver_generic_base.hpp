
#ifndef SVD_SOLVER_GENERIC_BASE_HPP_
#define SVD_SOLVER_GENERIC_BASE_HPP_

#include "svd_solver_traits.hpp"


namespace svd
{
    
template<typename derived_type>
class solverGenericBase
{
public:
  using sc_t = typename svd::details::traits<derived_type>::scalar_t;
  using wrap_solver_t = typename svd::details::traits<derived_type>::wrapped_solver_t;
  using der_t = derived_type;

  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };
  
  wrap_solver_t const & getConstRefToWrappedSolver() const {
    return this->underlying().getConstRefToWrappedSolverImpl();
  };
  
};
    
} // end namespace 

#endif
