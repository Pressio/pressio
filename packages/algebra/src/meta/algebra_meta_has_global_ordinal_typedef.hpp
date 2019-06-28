
#ifndef ALGEBRA_META_HAS_GLOBAL_ORDINAL_TYPEDEF_HPP_
#define ALGEBRA_META_HAS_GLOBAL_ORDINAL_TYPEDEF_HPP_

#include <type_traits>

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct has_global_ordinal_typedef : std::false_type{};

template <typename T>
struct has_global_ordinal_typedef<T,
		   typename
		    std::enable_if<
		     !std::is_void<typename T::global_ordinal_type
				   >::value
		     >::type
		   > : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
