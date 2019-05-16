
#ifndef CORE_META_HAS_UPDATE_OP_TYPEDEF_HPP_
#define CORE_META_HAS_UPDATE_OP_TYPEDEF_HPP_

#include <type_traits>

namespace rompp{ namespace core{ namespace meta {

/*
 * detect is a type T has a typedef named update_op
 */

template <typename T, typename = void>
struct has_update_op_typedef{
	static constexpr bool value = false;
	using type = void;
};

template <typename T>
struct has_update_op_typedef<
  T, mpl::enable_if_t<
       !std::is_void<
	 typename T::update_op
	 >::value
       >
  >{
	static constexpr bool value = true;
	using type = typename T::update_op;
};

}}} // namespace rompp::core::meta
#endif
