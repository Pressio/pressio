
#ifndef CORE_META_META_DETECTION_IDIOM_HPP 
#define CORE_META_META_DETECTION_IDIOM_HPP

#include <type_traits>

namespace core {
namespace meta {

// Source: cppreference. 
  
struct nonesuch {
  nonesuch() = delete;
  ~nonesuch() = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};


template <class Default,
	  class AlwaysVoid,
	  template <class...> class Op,
	  class... Args>
struct detector {
  using value_t = std::false_type;
  using type = Default;
};


template <typename... Ts>
struct my_make_void {
  using type = void;
};

template <typename... Ts>
using my_void_t = typename my_make_void<Ts...>::type;


template <class Default,
	  template <class...> class Op,
	  class... Args>
struct detector<Default, my_void_t<Op<Args...>>, Op, Args...> {
  using value_t = std::true_type;
  using type = Op<Args...>;
};


template <template <class...> class Op,
	  class... Args>
using is_detected =
  typename detector<nonesuch, void, Op, Args...>::value_t;


template <template <class...> class Op,
	  class... Args>
using detected_t =
  typename detector<nonesuch, void, Op, Args...>::type;


} // end namespace meta
} // end namespace core

#endif
