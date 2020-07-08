

#ifndef solvers_has_jacobian_typedef_HPP_
#define solvers_has_jacobian_typedef_HPP_

namespace pressio{ namespace solvers{ namespace predicates {

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

}}}
#endif
