
#ifndef ROM_MATRIX_PSEUDO_INVERSE_HPP_
#define ROM_MATRIX_PSEUDO_INVERSE_HPP_

#include "../rom_ConfigDefs.hpp"
#include "../../../CORE_MATRIX"
#include "../../../SVD_BASIC"

namespace rompp{
namespace rom{
namespace exp{


#ifdef HAVE_TRILINOS
template <typename mat_type,
	  typename std::enable_if<
	    core::details::traits<mat_type>::isEpetra &&
	    core::details::traits<mat_type>::isDense
	    >::type * = nullptr
	  >
auto pseudoInverse(const mat_type & A)
{
  const int nR = A.globalRows();
  const int nC = A.globalCols();
  Epetra_Map newMap(nC, 0, A.commCRef());
  mat_type C(newMap, nR);
  C.setZero();
  return C;
}
#endif 

  
}//end namespace exp
}//end namespace rom
}//end namespace rompp
#endif 
