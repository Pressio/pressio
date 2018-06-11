
#ifndef CORE_META_BASIC_HPP_
#define CORE_META_BASIC_HPP_

#include <type_traits>
#include <complex>

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

  template <typename... >
  using void_t = void;
  
  /////////////////////////////////////////////////

  // check if a type is default constructible
  // we leave this commented out, for now.
  // we use the std method instead.

  template<typename T>
  struct is_default_constructible
    : std::is_default_constructible<T> {};
  
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

  template <typename T,
	    typename enable = void>
  struct is_stdComplex : std::false_type{};

  template <typename T>
  struct is_stdComplex<T, typename
		       std::enable_if<
			    std::is_same<T,
					 std::complex<typename
						      T::value_type
						      >
					 >::value
			    >::type
		       > : std::true_type{};

  
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
