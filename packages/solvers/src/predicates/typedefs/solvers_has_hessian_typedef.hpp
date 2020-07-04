

#ifndef solvers_has_hessian_typedef_HPP_
#define solvers_has_hessian_typedef_HPP_

namespace pressio{ namespace solvers{ namespace meta {

template <typename T, typename enable = void>
struct has_hessian_typedef : std::false_type{};

template <typename T>
struct has_hessian_typedef<
  T,
  mpl::enable_if_t<
    !std::is_void<
      typename T::hessian_type
      >::value
    >
  > : std::true_type{};

}}}
#endif
