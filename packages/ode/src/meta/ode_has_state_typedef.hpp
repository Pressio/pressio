
#ifndef ODE_META_HAS_STATE_TYPEDEF_HPP_
#define ODE_META_HAS_STATE_TYPEDEF_HPP_

#include <type_traits>

namespace pressio{ namespace ode{ namespace meta {

template <typename T, typename enable = void>
struct has_state_typedef : std::false_type{};

template <typename T>
struct has_state_typedef<
  T,
  mpl::enable_if_t<
    !std::is_void<
      typename T::state_type
      >::value
    >
  > : std::true_type{};

}}}//end namespace pressio::ode::meta
#endif
