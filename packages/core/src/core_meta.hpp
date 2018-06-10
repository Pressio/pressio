
#ifndef CORE_META_HPP_
#define CORE_META_HPP_

#include <type_traits>


namespace core{
namespace meta {

  template<typename T>
  struct remove_const: std::remove_const<T>{};

  template<typename T>
  struct remove_reference : std::remove_reference<T>{};

  template<typename T>
  struct remove_pointer : std::remove_pointer<T>{};

  template<typename T>
  struct is_arithmetic : std::is_arithmetic<T>{};

  template<typename T>
  struct is_integral: std::is_integral<T>{};

  /////////////////////////////////////////////////

  // check if a type is default constructible
  // we leave this commented out, for now.
  // we use the std method instead.

  template<typename T>
  struct is_default_constructible : std::is_default_constructible<T> {};
  
  // template<typename T>
  // class is_default_constructible{
  //   typedef char yes;
  //   typedef struct { char arr[2]; } no;
  //   template<typename U>
  //   static decltype(U(), yes()) test(int);
  //   template<typename>
  //   static no test(...);
  // public:
  //   static const bool value = sizeof(test<T>(0)) == sizeof(yes);
  // };

  /////////////////////////////////////////////////

  // meta functions for checking if a
  // type has a public typedef named: scalar_type

  template <typename T, typename enable = void>
  struct has_scalarTypedef : std::false_type{};

  template <typename T>
  struct has_scalarTypedef<T,
			   typename std::enable_if<!std::is_void<typename T::scalar_type>::value
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

  // meta functions for checking if a
  // type has a public typedef named: ordinal_type

  template <typename T, typename enable = void>
  struct has_ordinalTypedef : std::false_type{};

  template <typename T>
  struct has_ordinalTypedef<T,
			   typename std::enable_if<!std::is_void<typename T::ordinal_type>::value
						   >::type
			   > : std::true_type{};
  
  /////////////////////////////////////////////////

  // meta functions for checking if a
  // type has a public typedef named: local_ordinal_type

  template <typename T, typename enable = void>
  struct has_localOrdinalTypedef : std::false_type{};

  template <typename T>
  struct has_localOrdinalTypedef<T,
			   typename std::enable_if<!std::is_void<typename T::local_ordinal_type>::value
						   >::type
			   > : std::true_type{};
  
  /////////////////////////////////////////////////

  // meta functions for checking if a
  // type has a public typedef named: global_ordinal_type

  template <typename T, typename enable = void>
  struct has_globalOrdinalTypedef : std::false_type{};

  template <typename T>
  struct has_globalOrdinalTypedef<T,
			   typename std::enable_if<!std::is_void<typename T::global_ordinal_type>::value
						   >::type
			   > : std::true_type{};
  
  /////////////////////////////////////////////////

  // meta functions for checking if a
  // type has a public typedef named: map_type

  template <typename T, typename enable = void>
  struct has_mapTypedef : std::false_type{};

  template <typename T>
  struct has_mapTypedef<T,
			   typename std::enable_if<!std::is_void<typename T::map_type>::value
						   >::type
			   > : std::true_type{};

  /////////////////////////////////////////////////

  // meta functions for checking if a
  // type has a public typedef named: comm_type

  template <typename T, typename enable = void>
  struct has_mpicommTypedef : std::false_type{};

  template <typename T>
  struct has_mpicommTypedef<T,
			    typename std::enable_if<!std::is_void<typename T::comm_type>::value
						   >::type
			   > : std::true_type{};
  
  /////////////////////////////////////////////////


} // namespace meta
} // namespace core
#endif





/////////////////////////////////////////////////
// template<class T1 , class T2> 
// struct same_instance_impl{ 
//   static bool same_instance( const T1& /* x1 */ , const T2& /* x2 */ ){
//     return false;
//   }
// };
// template< class T > 
// struct same_instance_impl<T,T>{ 
//   static bool same_instance( const T &x1 , const T &x2 ){
//     return (&x1 == &x2);
//   }
// };
// template< class T1 , class T2 > 
// bool same_instance( const T1 &x1 , const T2 &x2 ){
//   return same_instance_impl< T1,T2 >::same_instance( x1 , x2 );
// }

/////////////////////////////////////////////////
/////////////////////////////////////////////////
