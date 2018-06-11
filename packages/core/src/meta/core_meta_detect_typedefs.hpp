
#ifndef CORE_META_DETECT_TYPEDEFS_HPP_
#define CORE_META_DETECT_TYPEDEFS_HPP_

#include "core_meta_basic.hpp"

namespace core{
namespace meta {
  
  /////////////////////////////////////////////////
  // check type has a public typedef named: scalar_type
  /////////////////////////////////////////////////

  template <typename T, typename enable = void>
  struct has_scalarTypedef : std::false_type{};

  template <typename T>
  struct has_scalarTypedef<T,
			   typename
			   std::enable_if<
			     !std::is_void<typename T::scalar_type
					   >::value
			     >::type
			   > : std::true_type{};
  
  // An alternative way to do same thing above:
  //------------------------------------------
  // template <typename T> struct has_scTypedef {
  // private:
  //   template <typename T1> static typename T1::scalar test(int);
  //   template <typename> static void test(...);
  // public:
  //     static constexpr auto value = !std::is_void<decltype(test<T>(0))>::value;
  // };

  /////////////////////////////////////////////////
  // type has a public typedef named: ordinal_type
  /////////////////////////////////////////////////

  template <typename T, typename enable = void>
  struct has_ordinalTypedef : std::false_type{};

  template <typename T>
  struct has_ordinalTypedef<T,
			   typename
			    std::enable_if<
			     !std::is_void<typename T::ordinal_type
					   >::value
			     >::type
			   > : std::true_type{};
  
  /////////////////////////////////////////////////
  // type has a public typedef named: local_ordinal_type
  /////////////////////////////////////////////////

  template <typename T, typename enable = void>
  struct has_localOrdinalTypedef : std::false_type{};

  template <typename T>
  struct has_localOrdinalTypedef<T,typename
				 std::enable_if<
				     !std::is_void<typename
						   T::local_ordinal_type
						   >::value
				     >::type
				 > : std::true_type{};
  
  /////////////////////////////////////////////////
  // type has a public typedef named: global_ordinal_type
  /////////////////////////////////////////////////

  template <typename T, typename enable = void>
  struct has_globalOrdinalTypedef : std::false_type{};

  template <typename T>
  struct has_globalOrdinalTypedef<T,
				  typename
				  std::enable_if<
				    !std::is_void<typename
						  T::global_ordinal_type
						  >::value
				    >::type
				  > : std::true_type{};

  
  /////////////////////////////////////////////////
  // type has a public typedef named: map_type
  /////////////////////////////////////////////////

  template <typename T, typename enable = void>
  struct has_dataMapTypedef : std::false_type{};

  template <typename T>
  struct has_dataMapTypedef<T,
			typename
			std::enable_if<
			  !std::is_void<typename T::data_map_type
					>::value
			  >::type
			> : std::true_type{};

  
  /////////////////////////////////////////////////
  // type has a public typedef named: comm_type
  /////////////////////////////////////////////////

  template <typename T, typename enable = void>
  struct has_mpicommTypedef : std::false_type{};

  template <typename T>
  struct has_mpicommTypedef<T,
			    typename
			    std::enable_if<
			      !std::is_void<typename
					    T::communicator_type
					    >::value
			      >::type
			    > : std::true_type{};
  

} // namespace meta
} // namespace core
#endif
