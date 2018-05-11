
#ifndef COMMONBASE_CONFIGDEFS_HPP
#define COMMONBASE_CONFIGDEFS_HPP

#include "commonBase_config.h"


/// \brief Namespace for Tpetra implementation details.
/// \warning Do NOT rely on the contents of this namespace.
namespace commonBase{
namespace details{
namespace defaultTypes {

// namespace DefaultTypes {
//! Default value of Scalar template parameter.
typedef double scalar_t;

//! Default value of LocalOrdinal template parameter.
typedef int local_ordinal_t;

/// \typedef global_ordinal_type
typedef int global_ordinal_t;

} // namespace DefaultTypes
} // namespace Details
} // end of Tpetra namespace


namespace commonBase
{
  constexpr int _SUCCESS = 0;
  constexpr int _FAILURE = 1;

/*
  // #define EPETRA_CHK_ERR(a) { { int epetra_err = a; \
  //                             if ((epetra_err < 0 && Epetra_Object::GetTracebackMode() > 0) || \
  //                                 (epetra_err > 0 && Epetra_Object::GetTracebackMode() > 1)) { \
  //                     Epetra_Object::GetTracebackStream() << "Epetra ERROR " << epetra_err << ", " \
  //                          << __FILE__ << ", line " << __LINE__ << std::endl; }\
  //                     if (epetra_err != 0) return(epetra_err);  }\
  //                  }
*/

}

#endif
