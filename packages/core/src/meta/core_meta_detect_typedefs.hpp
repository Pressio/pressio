
#ifndef CORE_META_DETECT_TYPEDEFS_HPP_
#define CORE_META_DETECT_TYPEDEFS_HPP_

#include "core_meta_basic.hpp"

namespace core{
namespace meta {
  
  /////////////////////////////////////////////////
  // check type has a public typedef named: scalar_type
  /////////////////////////////////////////////////

  template <typename T, typename enable = void>
  struct has_scalar_typedef : std::false_type{};

  template <typename T>
  struct has_scalar_typedef<T,
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
  struct has_ordinal_typedef : std::false_type{};

  template <typename T>
  struct has_ordinal_typedef<T,
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
  struct has_local_ordinal_typedef : std::false_type{};

  template <typename T>
  struct has_local_ordinal_typedef<T,typename
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
  struct has_global_ordinal_typedef : std::false_type{};

  template <typename T>
  struct has_global_ordinal_typedef<T,
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
  struct has_data_map_typedef : std::false_type{};

  template <typename T>
  struct has_data_map_typedef<T,
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
  struct has_mpi_comm_typedef : std::false_type{};

  template <typename T>
  struct has_mpi_comm_typedef<T,
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
