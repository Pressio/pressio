
#ifndef CORE_META_HAS_ORDINAL_TYPEDEF_HPP_
#define CORE_META_HAS_ORDINAL_TYPEDEF_HPP_

#include <type_traits>


namespace rompp{ namespace core{ namespace meta {

  /////////////////////////////////////////////////
  // type has a public typedef named: ordinal_type
  /////////////////////////////////////////////////

  template <typename T, typename enable = void>
  struct has_ordinal_typedef : std::false_type{};

  template <typename T>
  struct has_ordinal_typedef<T,
			   typename
			    std::enable_if<
			     !std::is_void<typename T::ordinal_type
					   >::value
			     >::type
			   > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
