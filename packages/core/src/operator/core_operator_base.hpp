
#ifndef CORE_OPERATOR_BASE_HPP_
#define CORE_OPERATOR_BASE_HPP_

#include "../core_forward_declarations.hpp"

namespace core{

template <typename derived_t>
class OperatorBase
  : private core::details::CrtpBase<OperatorBase<derived_t>>{

public:
  // A_ (whatever that is) acts on X stored in Y
  template <typename operand_type,
	    typename result_type>
  void apply(const operand_type & X,
	     result_type & Y,
	     bool useTranspose = false)
  {
    this->underlying().applyImpl(X,Y,useTranspose);
  }

private:
  friend derived_t;
  friend core::details::CrtpBase<OperatorBase<derived_t>>;

  OperatorBase() = default;
  ~OperatorBase() = default;

};//end class

} // end namespace core
#endif
