
#ifndef CORE_META_HAS_SCALAR_TYPEDEF_HPP_
#define CORE_META_HAS_SCALAR_TYPEDEF_HPP_

#include <type_traits>

namespace rompp{ namespace core{ namespace meta {

  /////////////////////////////////////////////////
  // check type has a public typedef named: scalar_type
  /////////////////////////////////////////////////

  // template <typename T, typename enable = void>
  // struct has_scalar_typedef : std::false_type{};

  // template <typename T>
  // struct has_scalar_typedef<T,
  // 			   typename
  // 			   std::enable_if<
  // 			     !std::is_void<typename T::scalar_type
  // 					   >::value
  // 			     >::type
  // 			   > : std::true_type{};


  template <typename T>
  using has_scalar_typedef = typename T::scalar_type;


  // An alternative way to do same thing above:
  //------------------------------------------
  // template <typename T> struct has_scTypedef {
  // private:
  //   template <typename T1> static typename T1::scalar test(int);
  //   template <typename> static void test(...);
  // public:
  //     static constexpr auto value = !std::is_void<decltype(test<T>(0))>::value;
  // };

}}}//end namespace rompp::core::meta
#endif
