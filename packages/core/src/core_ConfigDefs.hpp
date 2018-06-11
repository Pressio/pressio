
#ifndef CORE_CONFIGDEFS_HPP_
#define CORE_CONFIGDEFS_HPP_

#include "core_config.h"
#include <stdlib.h>
#include <iostream>

namespace core{
namespace details {

  template<typename T, typename enable = void>
  struct traits;

  // traits<const T> == traits<T> 
  // For example, traits<Vector<const T> > != traits<Vector<T> >, but
  //              traits<const Vector<T> > == traits<Vector<T> >
  template<typename T> 
  struct traits<const T> : traits<T> {};

} // end namespace details


namespace defaultTypes {

  //! Default value of Scalar template parameter.
  using scalar_t = double;

  //! Default value of LocalOrdinal template parameter.
  using local_ordinal_t = int;

  /// default global_ordinal_type
  using global_ordinal_t = int;

  /// default type for error codes
  using errcode_t = int;

  // admissible types for epetra vector
  using epetra_scalar_t = double;
  using epetra_lo_t = int;
  using epetra_go_t1 = int;
  using epetra_go_t2 = long long;

} // namespace defaultTypes

constexpr defaultTypes::errcode_t _SUCCESS = 0;
constexpr defaultTypes::errcode_t _FAILURE = 1;
  
} // end of core namespace

#endif



/*
  // #define EPETRA_CHK_ERR(a) { { int epetra_err = a; \
  //                             if ((epetra_err < 0 && Epetra_Object::GetTracebackMode() > 0) || \
  //                                 (epetra_err > 0 && Epetra_Object::GetTracebackMode() > 1)) { \
  //                     Epetra_Object::GetTracebackStream() << "Epetra ERROR " << epetra_err << ", " \
  //                          << __FILE__ << ", line " << __LINE__ << std::endl; }\
  //                     if (epetra_err != 0) return(epetra_err);  }\
  //                  }
*/
