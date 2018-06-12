
#ifndef CORE_OPERATORS_BASE_HPP_
#define CORE_OPERATORS_BASE_HPP_

namespace core{

   
template<typename derived_type>
class arithmeticOperatorsBase{ 
public:
  derived_type operator+(const derived_type & other);
  derived_type operator-(const derived_type & other);
  derived_type operator*(const derived_type & other);
  
private:
  friend derived_type;
  arithmeticOperatorsBase() = default;
  ~arithmeticOperatorsBase() = default;

};//end class   

  

template<typename derived_type>
class compoundAssignmentOperatorsBase{
public:
  derived_type & operator+=(const derived_type & other);  
  derived_type & operator-=(const derived_type & other);

private:
  friend derived_type;
  compoundAssignmentOperatorsBase() = default;
  ~compoundAssignmentOperatorsBase() = default;
  
};//end class   


  
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



