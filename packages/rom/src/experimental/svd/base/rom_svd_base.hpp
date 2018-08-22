
#ifndef ROM_SVD_BASE_HPP_
#define ROM_SVD_BASE_HPP_

#include "rom_ConfigDefs.hpp"
#include "../rom_svd_traits.hpp"

namespace rom{

template<typename derived_type>
class SVDSolverBase
  : private core::details::CrtpBase<SVDSolverBase<derived_type>>
{
private:
  using traits = rom::details::traits<derived_type>;
  using mat_t = typename traits::matrix_t;

public:
  template<typename mat_t>
  void solve(mat_t & A){
    this->underlying().solveImpl(A);
  }

private:    
  SVDSolverBase() = default;
  ~SVDSolverBase() = default;

private:
  friend stepper_type;
  friend core::details::CrtpBase<SVDSolverBase<stepper_type>>;
  
};//end class

}//end namespace
#endif
