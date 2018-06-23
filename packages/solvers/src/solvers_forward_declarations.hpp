
#ifndef SOLVERS_FORWARD_DECLARATIONS_HPP_
#define SOLVERS_FORWARD_DECLARATIONS_HPP_

#include "solvers_ConfigDefs.hpp"

namespace solvers {






namespace experimental
{

// this is just a hack to have a solver to see 
// how overall workflow fits together
template<typename matrix_type, 
	 typename rhs_type, 
	 typename result_type,
	 // or an enun class, or similar
	 int algo_type = 0, 
	 typename Enable = void
	 >
class linearSolver;   

} // end experimental



} // end solvers
#endif
