
#ifndef CORE_META_DETECT_OPERATORS_HPP_
#define CORE_META_DETECT_OPERATORS_HPP_

#include "core_meta_basic.hpp"

namespace core{
namespace meta {

  /////////////////////////////////////////////////
  // detect [] operator
  /////////////////////////////////////////////////
  
  template<typename T,
	   typename ord_t,
	   typename enable = void>
  struct has_subscript_op : std::false_type { };

  template<typename T,
	   typename ord_t>
  struct has_subscript_op<T,
			  ord_t,
			  typename
			  std::enable_if<
			    !std::is_void<decltype(
				std::declval<T>()[std::declval<ord_t>()])
					  >::value, void
					>::type
			  > : std::true_type { };

  // void_t<
  // decltype(std::declval<T>()[std::declval<ord_t>()])
  //   >
  
  /////////////////////////////////////////////////
  // detect + operator
  /////////////////////////////////////////////////
  
  template<typename T,
	   typename U=T,
	   typename enable = void>
  struct has_add_op : std::false_type { };

  template<typename T,
	   typename U>
  struct has_add_op<T,U,
		    typename
		    std::enable_if<
		      !std::is_void<decltype(std::declval<T>() +
					     std::declval<U>())
				    >::value
		                  >::type
		    > : std::true_type{};

  
  /////////////////////////////////////////////////
  // detect - operator
  /////////////////////////////////////////////////
  
  template<typename T,
	   typename U=T,
	   typename enable = void>
  struct has_diff_op : std::false_type { };

  template<typename T, typename U>
  struct has_diff_op<T,U,
		    typename
		    std::enable_if<
		      !std::is_void<decltype(std::declval<T>() -
					     std::declval<U>())
				    >::value
		                  >::type
		    > : std::true_type{};

  
  /////////////////////////////////////////////////
  // detect * operator
  /////////////////////////////////////////////////
  
  template<typename T,
	   typename U=T,
	   typename enable = void>
  struct has_star_op : std::false_type { };

  template<typename T, typename U>
  struct has_star_op<T,U,
		    typename
		    std::enable_if<
		      !std::is_void<decltype(std::declval<T>() *
					     std::declval<U>())
				    >::value
		                  >::type
		    > : std::true_type{};


  /////////////////////////////////////////////////
  // detect += operator
  /////////////////////////////////////////////////
  
  template<typename T,
	   typename U=T,
	   typename enable = void>
  struct has_comp_assign_plus_op : std::false_type { };

  template<typename T, typename U>
  struct has_comp_assign_plus_op<T,U,
			    typename
			    std::enable_if<
			      !std::is_void<decltype(std::declval<T>() +=
						     std::declval<U>())
					    >::value
					  >::type
			    > : std::true_type{};


  /////////////////////////////////////////////////
  // detect -= operator
  /////////////////////////////////////////////////
  
  template<typename T,
	   typename U=T,
	   typename enable = void>
  struct has_comp_assign_minus_op : std::false_type { };

  template<typename T, typename U>
  struct has_comp_assign_minus_op<T,U,
			    typename
			    std::enable_if<
			      !std::is_void<decltype(std::declval<T>() -=
						     std::declval<U>())
					    >::value
					  >::type
			    > : std::true_type{};
  
 
} // namespace meta
} // namespace core
#endif
