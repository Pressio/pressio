
#ifndef CONTAINERS_META_HAS_DATA_MAP_TYPEDEF_HPP_
#define CONTAINERS_META_HAS_DATA_MAP_TYPEDEF_HPP_

#include <type_traits>

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct has_data_map_typedef : std::false_type{};

template <typename T>
struct has_data_map_typedef<T,
		typename
		std::enable_if<
		  !std::is_void<typename T::data_map_type
				>::value
		  >::type
		> : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
