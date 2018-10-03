
#ifndef ODE_CONFIGDEFS_HPP_
#define ODE_CONFIGDEFS_HPP_

#include "ode_config.h"
#include "../../core/src/core_ConfigDefs.hpp"
#include "../../core/src/meta/core_meta_basic.hpp"
#include "../../core/src/vector/core_vector_traits.hpp"
#include "../../core/src/matrix/core_matrix_traits.hpp"
#include "ode_enum_steppers.hpp"

namespace rompp{ namespace ode{

namespace details {
  template<typename T, typename enable = void>
  struct traits : core::details::traits<T> {};  
}// end namespace details

namespace impl {
  struct invalid_base_for_stepper{};
}// end namespace impl

} // end ode
}//end namespace rompp
#endif
