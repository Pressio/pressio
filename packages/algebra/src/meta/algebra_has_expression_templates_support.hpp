
#ifndef ALGEBRA_HAS_EXPRESSION_TEMPLATES_SUPPORT_HPP_
#define ALGEBRA_HAS_EXPRESSION_TEMPLATES_SUPPORT_HPP_

#include "../algebra_forward_declarations.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct has_expression_templates_support : std::false_type{};

template <typename T>
struct has_expression_templates_support
<T, ::rompp::mpl::enable_if_t<
      algebra::details::traits<T>::is_admissible_for_expression_templates
      >
 > : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
