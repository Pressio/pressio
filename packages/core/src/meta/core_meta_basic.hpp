
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

  //////////////////////////////////////////////////

  // template<template<class> class T, class U>
  // struct isDerivedFrom
  // {
  //   static constexpr bool value =
  //     decltype(isDerivedFrom::test(std::declval<U>()))::value;
  // private:
  //   template<class V>
  //   static decltype(static_cast<T<V>>(std::declval<U>()),
  // 		    std::true_type{}) test(const T<V>&);

  //   static std::false_type test(...);
  // };

  template<typename T, typename base_t>
  struct publiclyInheritsFrom : std::is_base_of<base_t,T>{};
  
  //////////////////////////////////////////////////

   
} // namespace meta
} // namespace core
#endif
