
#ifndef SOLVERS_META_BASIC_META_HPP_
#define SOLVERS_META_BASIC_META_HPP_

#include "../solvers_ConfigDefs.hpp"
#include "../../../mpl/src/detection_idiom.hpp"

namespace rompp{ namespace solvers{ namespace meta {

template <typename T>
using has_state_typedef = typename T::state_type;

template <typename T>
using has_residual_typedef = typename T::residual_type;

template <typename T>
using has_jacobian_typedef = typename T::jacobian_type;

template <typename T>
using has_scalar_typedef = typename T::scalar_type;

template <typename T>
using has_matrix_typedef = typename T::matrix_type;


}}} // namespace rompp::solvers::meta
#endif
