
#ifndef CORE_CONFIGDEFS_HPP_
#define CORE_CONFIGDEFS_HPP_

#include "core_crtp_helper.hpp"
#include "core_config.h"
#include <type_traits>
#include "core_shared_traits.hpp"
#include "meta/core_meta_basic.hpp"
#include "core_default_types.hpp"

namespace rompp{ namespace core{

namespace constants{

  // this is typically used as a template parameter
  // a positive quantity (e.g., a size of a vector)
  // is not known at compile-time, and the value
  // is defined at runtime
  constexpr int dynamic = -1;

}//end namespace constants


namespace details {

template<typename T, typename enable = void>
struct traits : public
containers_shared_traits<void, void,
			 false, false, false,
			 WrappedPackageIdentifier::Undefined,
			 false, false>{};

template<typename T>
struct traits<const T> : traits<T> {};

} // end namespace details
//--------------------------------------------

namespace exprtemplates{

struct plus_{
  template <typename a_t, typename b_t>
  auto operator()(const a_t & a, const b_t & b) const
  -> decltype(a+b) {
    return a + b;
  }
};

struct subtract_{
  template <typename a_t, typename b_t>
  auto operator()(const a_t & a, const b_t & b) const
  -> decltype(a-b) {
    return a - b;
  }
};

struct times_{
  template <typename a_t, typename b_t>
  auto operator()(const a_t & a, const b_t & b) const
  -> decltype(a*b) {
    return a * b;
  }
};

} // end namespace exprtemplates
//--------------------------------------------

}} // end of namespace rompp::core
#endif
