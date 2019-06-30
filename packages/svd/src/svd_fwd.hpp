
#ifndef SVD_FORWARD_DECLARATIONS_HPP_
#define SVD_FORWARD_DECLARATIONS_HPP_

#include "svd_ConfigDefs.hpp"

namespace rompp{ 
namespace svd {
    
template <typename matrix_type,
	  template<typename...> class lsv_type,
	  template<typename...> class rsv_type,
	  typename sval_type,
	  typename enable = void>
class Solver;  
  
} // end namespace svd
}//end namespace rompp
#endif
