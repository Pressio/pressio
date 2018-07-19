
#ifndef ROM_WEIGHTING_OPERATOR_BASE_HPP_
#define ROM_WEIGHTING_OPERATOR_BASE_HPP_

#include "rom_ConfigDefs.hpp"

namespace rom{
namespace experimental{

template <typename derived_t>
class weightingOperatorBase
  : private core::details::CrtpBase<weightingOperatorBase<derived_t>>
{

public:
  template <typename T>
  void leftMultiply(const T & yin, T & yout) 
  {
    this->underlying().leftMultiplyImpl(yin, yout);
  }

  template <typename T>
  void leftRightMultiply(const T & Min, T & Mout) 
  {
    this->underlying().leftRightMultiplyImpl(Min, Mout);
  }
  
private:
  friend derived_t;
  friend core::details::CrtpBase<weightingOperatorBase<derived_t>>;

  weightingOperatorBase() = default;
  ~weightingOperatorBase() = default;

};//end class

}//end namespace exp
}//end namespace rom
#endif
