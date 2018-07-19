
#ifndef CORE_OPERATORS_BASE_HPP_
#define CORE_OPERATORS_BASE_HPP_

#include <type_traits>
#include <iostream>

namespace core{
  
template<typename derived_type>
class ArithmeticOperatorsBase{ 
public:
  // disable (for now) if callers'type does not match input type
  template<typename U = derived_type,
  	   typename std::enable_if<
	     std::is_same<derived_type,U>::value
	     >::type * = nullptr>
  derived_type operator+(const U & other) const;
  //--------------------------------------------------
  // disable (for now) if callers'type does not match input type
  template<typename U = derived_type,
  	   typename std::enable_if<
	     std::is_same<derived_type,U>::value
	     >::type * = nullptr>
  derived_type operator-(const U & other) const;
  //--------------------------------------------------
  // disable (for now) if callers'type does not match input type
  template<typename U = derived_type,
  	   typename std::enable_if<
	     std::is_same<derived_type,U>::value
	     >::type * = nullptr>
  derived_type operator*(const U & other) const;
private:
  friend derived_type;

  ArithmeticOperatorsBase() = default;
  ~ArithmeticOperatorsBase() = default;

};//end class   

  
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
  

template<typename derived_type>
class CompoundAssignmentOperatorsBase{
public:
  // disable (for now) if callers'type does not match input type
  template<typename U = derived_type,
  	   typename std::enable_if<
	     std::is_same<derived_type,U>::value
	     >::type * = nullptr>
  derived_type & operator+=(const U & other);  
  //--------------------------------------------------
  // disable (for now) if callers'type does not match input type
  template<typename U = derived_type,
  	   typename std::enable_if<
	     std::is_same<derived_type,U>::value
	     >::type * = nullptr>
  derived_type & operator-=(const U & other);

private:
  friend derived_type;

  CompoundAssignmentOperatorsBase() = default;
  ~CompoundAssignmentOperatorsBase() = default;
  
};//end class   

  
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

  
template<typename derived_type,
	 typename scalar_type,
	 typename ordinal_type>
class SubscriptingOperatorsBase{
public:
  scalar_type & operator[] (ordinal_type index);
  scalar_type const & operator[] (ordinal_type index) const;

private:
  friend derived_type;

  SubscriptingOperatorsBase() = default;
  ~SubscriptingOperatorsBase() = default;  
 
};//end class   

} // end namespace core
#endif

