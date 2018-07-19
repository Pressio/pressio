
#ifndef ROM_SAMPLING_OPERATOR_COLLOCATION_HPP_
#define ROM_SAMPLING_OPERATOR_COLLOCATION_HPP_

#include "rom_ConfigDefs.hpp"
#include "./rom_sampling_operator_base.hpp"
#include <utility>

namespace rom{
namespace experimental{

template<typename matrix_type/*,
	  typename
	  std::enable_if<
	    core::meta::is_coreMatrixWrapper<matrix_type>::value && 
	    core::details::traits<matrix_type>::isEigen
	    >::type * = nullptr */
	  >
class samplingOperatorCollocation
  : public samplingOperatorBase<samplingOperatorCollocation<matrix_type>>
{

public:

  template <typename ... Args>
  samplingOperatorCollocation(Args&&... rest)
    : P_( std::forward<Args>(rest)... ){}

  ~samplingOperatorCollocation() = default;
  
private:
  template <typename vector_type>
  void applyImpl(const vector_type & yin,
		 vector_type & yout) const{
    //    *yout.data()
  }

private:
  matrix_type P_;

private:
  friend samplingOperatorBase<samplingOperatorCollocation<matrix_type>>;

};//end class

}//end namespace exp
}//end namespace rom
#endif
