
#ifndef CORE_META_META_DETECTION_IDIOM_HPP 
#define CORE_META_META_DETECTION_IDIOM_HPP

#include "core_meta_basic.hpp"

namespace rompp{
namespace core {
namespace meta {

// Source: cppreference. 
  
struct nonesuch {
  nonesuch() = delete;
  ~nonesuch() = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};

// already defined in core_meta_basic
// template <typename...>
// using my_void_t = void;

template <class Default,
	  class AlwaysVoid,
	  template <class...> class Op,
	  class...>
struct detector {
  using value_t = std::false_type;
  using type = Default;
};

template <class Default,
	  template <class...> class Op,
	  class... Args>
struct detector<Default, void_t<Op<Args...>>, Op, Args...> {
  using value_t = std::true_type;
  using type = Op<Args...>;
};

//================================================
  
template <template <class...> class Op, class... Args>
using is_detected = typename detector<nonesuch, void, Op, Args...>::value_t;

template <template <class...> class Op, class... Args>
using detected_t = typename detector<nonesuch, void, Op, Args...>::type;
  
template <class T, template<class...> class Op, class... Args>
using detected_or = detector<T, void, Op, Args...>;  

template <class T, template<class...> class Op, class... Args>
using detected_or_t = typename detected_or<T, Op, Args...>::type;

template <class T, template<class...> class Op, class... Args>
using is_detected_exact = std::is_same<T, detected_t<Op, Args...>>;
  
} // end namespace meta
} // end namespace core

}//end namespace rompp
#endif
