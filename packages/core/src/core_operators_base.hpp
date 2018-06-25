
#ifndef CORE_OPERATORS_BASE_HPP_
#define CORE_OPERATORS_BASE_HPP_

#include <type_traits>
#include <iostream>

namespace core{
  
template<typename derived_type>
class arithmeticOperatorsBase{ 
public:
  // disable (for now) if callers'type does not match input type
  template<typename U = derived_type,
  	   typename std::enable_if<
	     std::is_same<derived_type,U>::value
	     >::type * = nullptr>
  derived_type operator+(const U & other);
  //--------------------------------------------------
  // disable (for now) if callers'type does not match input type
  template<typename U = derived_type,
  	   typename std::enable_if<
	     std::is_same<derived_type,U>::value
	     >::type * = nullptr>
  derived_type operator-(const U & other);
  //--------------------------------------------------
  // disable (for now) if callers'type does not match input type
  template<typename U = derived_type,
  	   typename std::enable_if<
	     std::is_same<derived_type,U>::value
	     >::type * = nullptr>
  derived_type operator*(const U & other);
private:
  friend derived_type;
  arithmeticOperatorsBase() = default;
  ~arithmeticOperatorsBase() = default;
};//end class   

  
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
  

template<typename derived_type>
class compoundAssignmentOperatorsBase{
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
  compoundAssignmentOperatorsBase() = default;
  ~compoundAssignmentOperatorsBase() = default;
  
};//end class   

  
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

  
template<typename derived_type,
	 typename scalar_type,
	 typename ordinal_type>
class subscriptingOperatorsBase{
public:
  scalar_type & operator[] (ordinal_type index);
  scalar_type const & operator[] (ordinal_type index) const;

private:
  friend derived_type;
  subscriptingOperatorsBase() = default;
  ~subscriptingOperatorsBase() = default;  
 
};//end class   


} // end namespace core
#endif



