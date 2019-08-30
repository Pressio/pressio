
#ifndef SOLVERS_IS_LEGITIMATE_LINE_SEARCH_TAG_HPP_
#define SOLVERS_IS_LEGITIMATE_LINE_SEARCH_TAG_HPP_

#include "../solvers_line_search_tags.hpp"

namespace pressio{ namespace solvers{ namespace meta {

template <typename T, typename enable = void>
struct is_legitimate_line_search_tag
  : std::false_type{};

template <typename T>
struct is_legitimate_line_search_tag<
  T,
  ::pressio::mpl::enable_if_t<
    std::is_same<
      T,
      ::pressio::solvers::iterative::gn::ArmijoLineSearch>::value
    >
  > : std::true_type{};

}}} // namespace pressio::solvers::meta
#endif
