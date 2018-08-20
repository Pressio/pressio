
#ifndef ROM_POD_HPP_
#define ROM_POD_HPP_

#include "rom_ConfigDefs.hpp"
#include "CORE_MATRIX"

namespace rom{
namespace exp{

template <typename mat_type,
	  typename std::enable_if<
	    core::details::traits<mat_type>::isEpetra &&
	    core::details::traits<mat_type>::isDense
	    >::type * = nullptr
	  >
auto pod(const mat_type & A)
{
  // const int nR = A.globalRows();
  // const int nC = A.globalCols();
  // // Epetra_Map newMap(nC, 0, A.commCRef());
  // // mat_type C(newMap, nR);
  // // C.setZero();
  // return C;

  // if sz is not None and engy is not None:
  //     print('Warning: Specified both sz and engy; sz ignored.')
  // if svd is None:
  //     from scipy.linalg import svd as svdloc
  //     svd = lambda Y: svdloc(Y, full_matrices=False, compute_uv=True)
  // U, s, V = svd(X)
  // rnk = sum(s>=rnkdef_tol)
  // if szbnds is None:
  //     szbnds = (1, rnk)
  // szbnds = (max(1, szbnds[0]), min(rnk, szbnds[1]))
  // if engy is not None:
  //     sz = rnk if engy >= 1.0 else np.nonzero(trunc_engy(s)>=engy)[0][0]+1
  // if sz is not None:
  //     sz = min(max(sz, szbnds[0]), szbnds[1])
  // return U[..., :sz], s[:sz], V[..., :sz]

  
  double a;
  return a;
}
//------------------------------------------------------

  
}//end namespace exp
}//end namespace rom
#endif 
