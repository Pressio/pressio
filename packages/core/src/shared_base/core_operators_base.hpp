
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
  
// template<typename derived_type,
// 	 typename scalar_type,
// 	 typename ordinal_type>
// class Subscripting1DOperatorsBase{
// public:
//   scalar_type & operator[] (ordinal_type index){
//     return static_cast<derived_type &>(*this)[index];
//   }

//   scalar_type const & operator[] (ordinal_type index) const{
//     return static_cast<const derived_type &>(*this)[index];
//   }

// private:
//   friend derived_type;
//   Subscripting1DOperatorsBase() = default;
//   ~Subscripting1DOperatorsBase() = default;  
 
// };//end class   
// //--------------------------------------------------


// template<typename derived_type,
// 	 typename scalar_type,
// 	 typename row_ordinal_type,
// 	 typename col_ordinal_type = row_ordinal_type>
// class Subscripting2DOperatorsBase{
// public:
//   scalar_type & operator()(row_ordinal_type irow, 
//                            col_ordinal_type icol);
    
//   scalar_type const & operator()(row_ordinal_type irow, 
//                                  col_ordinal_type icol) const;

// private:
//   friend derived_type;
//   Subscripting2DOperatorsBase() = default;
//   ~Subscripting2DOperatorsBase() = default;  
 
// };//end class   
// //--------------------------------------------------


} // end namespace core
}//end namespace rompp
#endif

