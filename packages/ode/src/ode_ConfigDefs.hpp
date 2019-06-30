
#ifndef ODE_CONFIGDEFS_HPP_
#define ODE_CONFIGDEFS_HPP_

#include "ode_config.h"
#include "../../containers/src/containers_ConfigDefs.hpp"
#include "../../containers/src/meta/containers_meta_basic.hpp"
#include "../../containers/src/vector/containers_vector_traits.hpp"
#include "../../containers/src/matrix/containers_matrix_traits.hpp"
#include "../../containers/src/multi_vector/containers_multi_vector_traits.hpp"
#include "../../containers/src/vector/containers_vector_meta.hpp"
#include "../../containers/src/multi_vector/containers_multi_vector_meta.hpp"
#include "../../containers/src/matrix/containers_matrix_meta.hpp"
#include "../../containers/src/meta/containers_is_wrapper.hpp"

#ifdef DEBUG_PRINT
#include "../../utils/src/io/utils_print_helper.hpp"
#endif

#include "ode_enum.hpp"

namespace rompp{ namespace ode{ namespace details {

template<typename T, typename enable = void> struct traits{};
template<typename T>  struct traits<const T> : traits<T> {};

}}}// end namespace rompp::ode::details
#endif
