
#ifndef ODE_META_HAS_JACOBIAN_TYPEDEF_HPP_
#define ODE_META_HAS_JACOBIAN_TYPEDEF_HPP_

#include <type_traits>

namespace pressio{ namespace ode{ namespace meta {

template <typename T, typename enable = void>
struct has_jacobian_typedef : std::false_type{};

template <typename T>
struct has_jacobian_typedef<
  T,
  mpl::enable_if_t<
    !std::is_void<
      typename T::jacobian_type
      >::value
    >
  > : std::true_type{};

}}}//end namespace pressio::ode::meta
#endif
