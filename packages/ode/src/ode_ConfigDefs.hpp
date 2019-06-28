
#ifndef ODE_CONFIGDEFS_HPP_
#define ODE_CONFIGDEFS_HPP_

#include "ode_config.h"
#include "../../algebra/src/algebra_ConfigDefs.hpp"
#include "../../algebra/src/meta/algebra_meta_basic.hpp"
#include "../../algebra/src/vector/algebra_vector_traits.hpp"
#include "../../algebra/src/matrix/algebra_matrix_traits.hpp"
#include "../../algebra/src/multi_vector/algebra_multi_vector_traits.hpp"
#include "../../algebra/src/vector/algebra_vector_meta.hpp"
#include "../../algebra/src/multi_vector/algebra_multi_vector_meta.hpp"
#include "../../algebra/src/matrix/algebra_matrix_meta.hpp"
#include "../../algebra/src/algebra_is_wrapper.hpp"

#ifdef DEBUG_PRINT
#include "../../utils/src/io/utils_print_helper.hpp"
#endif

#include "ode_enum.hpp"

namespace rompp{ namespace ode{ namespace details {

template<typename T, typename enable = void> struct traits{};
template<typename T>  struct traits<const T> : traits<T> {};

}}}// end namespace rompp::ode::details
#endif
