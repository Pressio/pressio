
#ifndef SVD_CONFIGDEFS_HPP
#define SVD_CONFIGDEFS_HPP

#include "svd_config.h"
#include "core_ConfigDefs.hpp"

namespace svd{
  
// put here definitions
enum class svdType{ 
  // /*FULL SVD
  //   The matrix U'n is thus m×m, Σn is m×n diagonal, and V is n×n.
  // */
  // full, 

  // /*THIN SVD: Only the n column vectors of U corresponding to the 
  //   row vectors of V* are calculated. The remaining column 
  //   vectors of U are not calculated. This is significantly 
  //   quicker and more economical than the full SVD if n << m. 
  //   The matrix U'n is thus m×n, Σn is n×n diagonal, and V is n×n.
  // */
  // thin,

  // /*Only the r column vectors of U and r row vectors of V* 
  //   corresponding to the non-zero singular values Σr are calculated. 
  //   The remaining vectors of U and V* are not calculated. 
  //   This is quicker and more economical than the thin SVD if r << n. 
  //   The matrix Ur is thus m×r, Σr is r×r diagonal, and Vr* is r×n.*/
  // compact,


  /*Only the t column vectors of U and t row vectors of V*
    corresponding to the t largest singular values Σt are calculated. 
    The rest of the matrix is discarded. This can be much quicker 
    and more economical than the compact SVD if t≪r. The matrix Ut 
    is thus m×t, Σt is t×t diagonal, and Vt* is t×n.*/
  truncated
  
};


  
//namespace details {  
// template<typename T, typename enable = void>
// struct traits : core::details::traits<T> {};
//} // end namespace details


} // end of svd namespace

#endif
