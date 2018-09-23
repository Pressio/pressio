
#ifndef ROM_INCREMENTAL_SOLUTION_BASE_HPP_
#define ROM_INCREMENTAL_SOLUTION_BASE_HPP_

#include "../rom_ConfigDefs.hpp"

namespace rompp{
namespace rom{
namespace exp{

template <template <typename...> class derived_type,
	  typename state_type,
	  typename ... Args>
class incrementalSolutionBase
{
public:
  static constexpr bool advanceIncrement = true;
  static constexpr bool advanceFull = !advanceIncrement;

protected:
  state_type const * y0ptr_;
  state_type yFull_;

private:
  using derived_t = derived_type<state_type, Args...>;
  friend derived_t; 

  incrementalSolutionBase(const state_type & y0)
    : y0ptr_(&y0), yFull_(y0){}

  ~incrementalSolutionBase() = default;

};//end class

}//end namespace polices
}//end namespace ode  
}//end namespace rompp
#endif 
