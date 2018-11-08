
#ifndef CORE_DEFAULT_TYPES_HPP_
#define CORE_DEFAULT_TYPES_HPP_

#include "core_config.h"

namespace rompp{ namespace core{ namespace default_types {

  //! Default value of Scalar template parameter.
  using scalar_t = double;

  //! Default value of LocalOrdinal template parameter.
  using local_ordinal_t = int;

  /// default global_ordinal_type
  using global_ordinal_t = int;

  // Unsigned int type
  using uint = unsigned int;

  // /// default type for error codes
  // using errcode_t = int;

#ifdef HAVE_TRILINOS
  // admissible types for epetra vector
  using epetra_scalar_t = double;
  using epetra_lo_t = int;
  using epetra_go_t1 = int;
  using epetra_go_t2 = long long;
#endif

    
}}} // end of namespace rompp::core::default_types
#endif
