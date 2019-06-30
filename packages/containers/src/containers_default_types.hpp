
#ifndef CONTAINERS_DEFAULT_TYPES_HPP_
#define CONTAINERS_DEFAULT_TYPES_HPP_

#include "containers_config.h"
#include "../../utils/src/utils_ConfigDefs.hpp"

namespace rompp{ namespace containers{ namespace default_types {

  //! Default value of Scalar template parameter.
  using scalar_t = double;

  //! Default value of LocalOrdinal template parameter.
  using local_ordinal_t = int;

  /// default global_ordinal_type
  using global_ordinal_t = int;

  // Unsigned int type
  using uint = unsigned int;

#ifdef HAVE_TRILINOS
  // admissible types for epetra vector
  using epetra_scalar_t = double;
  using epetra_lo_t = int;
  using epetra_go_t1 = int;
  using epetra_go_t2 = long long;
#endif


}}} // end of namespace rompp::containers::default_types
#endif
