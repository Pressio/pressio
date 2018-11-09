
#ifndef ROM_OPERATOR_BASE_HPP_
#define ROM_OPERATOR_BASE_HPP_

#include "../rom_ConfigDefs.hpp"

namespace rompp{ namespace rom{

template <typename derived_t>
class OperatorBase
  : private core::details::CrtpBase<OperatorBase<derived_t>>{

public:

  // A (whatever that is) acts on X, return result 
  template <typename operand_type>
  auto apply(const operand_type & X)
  -> decltype(this->underlying().applyImpl(X))
  {
    return this->underlying().applyImpl(X);
  }

  // A (whatever that is) acts on X
  // result stored in Y
  template <typename operand_t1,
	    typename operand_t2>
  void apply(const operand_t1 & X,
	     operand_t2 & Y){
    this->underlying().applyImpl(X,Y);
  }

  // A (whatever that is) acts on X from right, return result
  // => X A 
  template <typename operand_type>
  auto applyRight(const operand_type & X)
  -> decltype(this->underlying().applyRightImpl(X))
  {
    return this->underlying().applyRightImpl(X);
  }
  
  // A^T (whatever that is) acts on X, return result
  template <typename operand_type>
  auto applyTranspose(const operand_type & X)
  -> decltype(this->underlying().applyTransposeImpl(X)){
    return this->underlying().applyTransposeImpl(X);
  }

  // A^T (whatever that is) acts on X,
  // result stored in Y
  template <typename operand_t1,
	    typename operand_t2>
  void applyTranspose(const operand_t1 & X,
		      operand_t2 & Y){
    this->underlying().applyTransposeImpl(X,Y);
  }
  
private:
  friend derived_t;
  friend core::details::CrtpBase<OperatorBase<derived_t>>;
  OperatorBase() = default;
  ~OperatorBase() = default;

};//end class

}} // end namespace rompp::rom
#endif
