
#ifndef PRESSIO_MPL_DETECTION_IDIOM_HPP 
#define PRESSIO_MPL_DETECTION_IDIOM_HPP

#include <type_traits>
#include "void_t.hpp"

namespace pressio{ namespace mpl{
  
struct nonesuch {
  nonesuch() = delete;
  ~nonesuch() = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};

// already defined in containers_meta_basic
// template <typename...>
// using my_void_t = void;

template <class Default,
	  class AlwaysVoid,
	  template <class...> class Op,
	  class...>
struct detector {
  constexpr static auto value = false;  
  using type = Default;
};

template <class Default,
	  template <class...> class Op,
	  class... Args>
struct detector<Default, 
                ::pressio::mpl::void_t<Op<Args...>>, Op, Args...> {
  constexpr static auto value = true;
  using type = Op<Args...>;
};

//================================================
  
template <template <class...> class Op, class... Args>
using is_detected = detector<nonesuch, void, Op, Args...>;

template <template <class...> class Op, class... Args>
using detected_t = typename detector<nonesuch, void, Op, Args...>::type;
  
template <class T, template<class...> class Op, class... Args>
using detected_or = detector<T, void, Op, Args...>;  

template <class T, template<class...> class Op, class... Args>
using detected_or_t = typename detected_or<T, Op, Args...>::type;

template <class T, template<class...> class Op, class... Args>
using is_detected_exact = std::is_same<T, detected_t<Op, Args...>>;
  

}}//end namespace pressio::mpl
#endif
