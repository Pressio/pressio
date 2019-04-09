
#ifndef SOLVERS_IS_NON_DEFAULT_LINE_SEARCH_TAG_HPP_
#define SOLVERS_IS_NON_DEFAULT_LINE_SEARCH_TAG_HPP_

#include "../solvers_line_search_tags.hpp"

namespace rompp{ namespace solvers{ namespace meta {

template <typename T, typename enable = void>
struct is_non_default_line_search_tag
  : std::false_type{};

template <typename T>
struct is_non_default_line_search_tag<
  T,
  core::meta::enable_if_t<
    std::is_same<
      T,
      ::rompp::solvers::iterative::gn::ArmijoLineSearch>::value
    >
  > : std::true_type{};

}}} // namespace rompp::solvers::meta
#endif
