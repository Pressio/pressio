
#ifndef ODE_CONFIGDEFS_HPP_
#define ODE_CONFIGDEFS_HPP_

#include "ode_config.h"
#include "../../core/src/core_ConfigDefs.hpp"
#include "../../core/src/meta/core_meta_basic.hpp"
#include "../../core/src/vector/core_vector_traits.hpp"
#include "../../core/src/matrix/core_matrix_traits.hpp"
#include "../../core/src/multi_vector/core_multi_vector_traits.hpp"
#include "../../core/src/vector/core_vector_meta.hpp"
#include "../../core/src/multi_vector/core_multi_vector_meta.hpp"
#include "../../core/src/matrix/core_matrix_meta.hpp"
#ifdef DEBUG_PRINT
#include "../../core/src/io/core_print_helper.hpp"
#endif

#include "ode_enum.hpp"

namespace rompp{ namespace ode{ namespace details {

template<typename T, typename enable = void> struct traits{};
template<typename T>  struct traits<const T> : traits<T> {};

}}}// end namespace rompp::ode::details
#endif
