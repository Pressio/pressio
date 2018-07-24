
#ifndef CORE_OPERATORS_BASE_HPP_
#define CORE_OPERATORS_BASE_HPP_

#include <type_traits>

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
class Subscripting1DOperatorsBase{
public:
  scalar_type & operator[] (ordinal_type index);
  scalar_type const & operator[] (ordinal_type index) const;

private:
  friend derived_type;
  Subscripting1DOperatorsBase() = default;
  ~Subscripting1DOperatorsBase() = default;  
 
};//end class   

  
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////


template<typename derived_type,
   typename scalar_type,
   typename ordinal_type>
class Subscripting2DOperatorsBase{
public:
  scalar_type & operator()(ordinal_type irow, 
                           ordinal_type icol);
  scalar_type const & operator()(ordinal_type irow, 
                                 ordinal_type icol) const;

private:
  friend derived_type;
  Subscripting2DOperatorsBase() = default;
  ~Subscripting2DOperatorsBase() = default;  
 
};//end class   


} // end namespace core
#endif

