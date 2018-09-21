
#ifndef CORE_CRTP_HELPER_HPP_
#define CORE_CRTP_HELPER_HPP_

#include <type_traits>

namespace core{
namespace details {  

//---------------------------------------
// CRTP HELPER BASE CLASS
//---------------------------------------
template <typename T, typename enable = void>
struct CrtpBase;
  
template <typename T,
	  template<typename, typename...> class crtpType,
	  typename ... Args>
struct CrtpBase< crtpType<T, Args...>>{
  T & underlying() {
    return static_cast<T&>(*this);
  }
  T const & underlying() const {
    return static_cast<T const&>(*this);
  }

private:
  CrtpBase() = default;
  ~CrtpBase() = default;
  friend crtpType<T, Args...>;
};//end class


template <typename T, int a, int b,
	  template<typename, int, int> class crtpType>
struct CrtpBase< crtpType<T, a, b> >{
  T & underlying() {
    return static_cast<T&>(*this);
  }
  T const & underlying() const {
    return static_cast<T const&>(*this);
  }

private:
  CrtpBase() = default;
  ~CrtpBase() = default;
  friend crtpType<T, a, b>;
};//end class

  
template <typename T, int a,
	  template<typename, int> class crtpType>
struct CrtpBase< crtpType<T, a> >{
  T & underlying() {
    return static_cast<T&>(*this);
  }
  T const & underlying() const {
    return static_cast<T const&>(*this);
  }

private:
  CrtpBase() = default;
  ~CrtpBase() = default;
  friend crtpType<T, a>;
};//end class

  
} // end namespace details
} // end of core namespace

#endif
