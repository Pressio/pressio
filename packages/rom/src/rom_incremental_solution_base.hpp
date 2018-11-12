
#ifndef ROM_INCREMENTAL_SOLUTION_BASE_HPP_
#define ROM_INCREMENTAL_SOLUTION_BASE_HPP_

#include "rom_ConfigDefs.hpp"

namespace rompp{ namespace rom{ namespace exp{

template <typename derived_type, typename state_type>
class IncrementalSolutionBase{
private:
  const state_type * y0FOM_;
  mutable state_type yFOM_;

private:
  friend derived_type; 
  IncrementalSolutionBase(const state_type & y0fom)
    : y0FOM_(&y0fom), yFOM_(y0fom){}
  ~IncrementalSolutionBase() = default;

};//end class

}}}//end namespace rompp::ode::polices
#endif 
