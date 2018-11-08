
#ifndef CORE_CONFIGDEFS_HPP_
#define CORE_CONFIGDEFS_HPP_

#include "core_crtp_helper.hpp"
#include "core_config.h"
#include <type_traits>
#include "core_shared_traits.hpp"
#include "meta/core_meta_basic.hpp"

namespace rompp{ namespace core{



namespace details {
  
template<typename T, typename enable = void>
struct traits : public
containers_shared_traits<void, void,
			 false, false, false,
			 WrappedPackageIdentifier::Undefined,
			 false>{};
  
template<typename T> 
struct traits<const T> : traits<T> {};
  
} // end namespace details
//--------------------------------------------

namespace exprtemplates{

struct plus_{
  template <typename a_t, typename b_t>
  auto operator()(const a_t & a, const b_t & b) const {
    return a + b;
  }  
};
  
struct subtract_{
  template <typename a_t, typename b_t>
  auto operator()(const a_t & a, const b_t & b) const {
    return a - b;
  }  
};

struct times_{
  template <typename a_t, typename b_t>
  auto operator()(const a_t & a, const b_t & b) const {
    return a * b;
  }  
};

} // end namespace exprtemplates
//--------------------------------------------
    
   
namespace default_types {

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

} // namespace default_types

    
}} // end of namespace rompp::core
#endif
