
#ifndef ROM_SAMPLING_OPERATOR_BASE_HPP_
#define ROM_SAMPLING_OPERATOR_BASE_HPP_

#include "rom_ConfigDefs.hpp"

namespace rom{
namespace experimental{

template <typename derived_t>
class samplingOperatorBase
  : private core::details::CrtpBase<samplingOperatorBase<derived_t>>
{

public:
  template <typename vector_type>
  void apply(const vector_type & yin,
	     vector_type & yout) const
  {
    this->underlying().applyImpl(yin, yout);
  }
  
private:
  friend derived_t;
  friend core::details::CrtpBase<samplingOperatorBase<derived_t>>;

  samplingOperatorBase() = default;
  ~samplingOperatorBase() = default;  

};//end class

}//end namespace exp
}//end namespace rom
#endif
