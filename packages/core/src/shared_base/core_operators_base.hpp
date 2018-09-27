
#ifndef CORE_OPERATORS_BASE_HPP_
#define CORE_OPERATORS_BASE_HPP_

#include <type_traits>

namespace rompp{
namespace core{
  
template<typename derived_type>
class ArithmeticOperatorsBase{ 
public:
  // disable (for now) if callers'type does not match input type
  template<typename T,
  	   typename std::enable_if<
	     std::is_same<derived_type,T>::value
	     >::type * = nullptr>
  derived_type operator+(const T & other) const;

  // disable (for now) if callers'type does not match input type
  template<typename T,
  	   typename std::enable_if<
	     std::is_same<derived_type,T>::value
	     >::type * = nullptr>
  derived_type operator-(const T & other) const;

  // disable (for now) if callers'type does not match input type
  template<typename T,
  	   typename std::enable_if<
	     std::is_same<derived_type,T>::value
	     >::type * = nullptr>
  derived_type operator*(const T & other) const;

private:
  friend derived_type;
  ArithmeticOperatorsBase() = default;
  ~ArithmeticOperatorsBase() = default;

};//end class
//--------------------------------------------------
  
} // end namespace core
}//end namespace rompp
#endif

