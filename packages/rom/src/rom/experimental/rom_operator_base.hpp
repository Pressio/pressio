
#ifndef ROM_OPERATOR_BASE_HPP_
#define ROM_OPERATOR_BASE_HPP_

#include "rom_ConfigDefs.hpp"

namespace rom{

template <typename derived_t>
class OperatorBase
  : private core::details::CrtpBase<OperatorBase<derived_t>>{

public:
  // // A_ (whatever that is) acts on X stored in Y
  // template <typename operand_type,
  // 	    typename result_type,
  // 	    bool transpose = false>
  // void apply(const operand_type & X,
  // 	     result_type & Y)
  // {
  //   this->underlying().template applyImpl<operand_type,
  // 					  result_type,
  // 					  transpose>(X,Y);
  // }

  // A_ (whatever that is) acts on X stored in Y
  template <typename operand_type>
  auto apply(const operand_type & X)
  {
    return this->underlying().applyImpl(X);
  }

  // A_ (whatever that is) acts on X stored in Y
  template <typename operand_type>
  auto applyTransp(const operand_type & X)
  {
    return this->underlying().applyTranspImpl(X);
  }
  
  
private:
  friend derived_t;
  friend core::details::CrtpBase<OperatorBase<derived_t>>;

  OperatorBase() = default;
  ~OperatorBase() = default;

};//end class

} // end namespace rom
#endif
