
#ifndef CONTAINERS_HAS_EXPRESSION_TEMPLATES_SUPPORT_HPP_
#define CONTAINERS_HAS_EXPRESSION_TEMPLATES_SUPPORT_HPP_

#include "../containers_fwd.hpp"

namespace rompp{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct has_expression_templates_support : std::false_type{};

template <typename T>
struct has_expression_templates_support
<T, ::rompp::mpl::enable_if_t<
      containers::details::traits<T>::is_admissible_for_expression_templates
      >
 > : std::true_type{};

}}}//end namespace rompp::containers::meta
#endif
