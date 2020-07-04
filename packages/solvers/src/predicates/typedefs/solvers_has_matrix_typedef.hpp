

#ifndef solvers_has_matrix_typedef_HPP_
#define solvers_has_matrix_typedef_HPP_

namespace pressio{ namespace solvers{ namespace meta {

template <typename T, typename enable = void>
struct has_matrix_typedef : std::false_type{};

template <typename T>
struct has_matrix_typedef<
  T,
  mpl::enable_if_t<
    !std::is_void<
      typename T::matrix_type
      >::value
    >
  > : std::true_type{};

}}}
#endif
