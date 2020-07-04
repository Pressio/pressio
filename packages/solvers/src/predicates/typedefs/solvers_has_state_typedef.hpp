

#ifndef solvers_has_state_typedef_HPP_
#define solvers_has_state_typedef_HPP_

namespace pressio{ namespace solvers{ namespace meta {

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

}}}
#endif
