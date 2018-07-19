
#ifndef ROM_BASIS_OPERATOR_BASE_HPP_
#define ROM_BASIS_OPERATOR_BASE_HPP_

#include "rom_ConfigDefs.hpp"

namespace rom{
namespace experimental{

template <typename derived_t>
class basisOperatorBase
  : private core::details::CrtpBase<basisOperatorBase<derived_t>>
{

public:
  template <typename T>
  void project(const T & yin, T & yout) const
  {
    this->underlying().projectImpl(yin, yout);
  }

  template <typename T>
  void leftMultiply(const T & yin, T & yout) 
  {
    this->underlying().leftMultiplyImpl(yin, yout);
  }

  template <typename T, typename T1>
  void rightMultiply(const T & yin, T1 & yout) 
  {
    this->underlying().rightMultiplyImpl(yin, yout);
  }
  
private:
  friend derived_t;
  friend core::details::CrtpBase<basisOperatorBase<derived_t>>;

  basisOperatorBase() = default;
  ~basisOperatorBase() = default;

};//end class

}//end namespace exp
}//end namespace rom
#endif
