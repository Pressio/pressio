
#ifndef ROM_SAMPLING_OPERATOR_IDENTITY_HPP_
#define ROM_SAMPLING_OPERATOR_IDENTITY_HPP_

#include "rom_ConfigDefs.hpp"
#include "./rom_sampling_operator_base.hpp"
#include <type_traits>

namespace rom{
namespace experimental{

template<typename matrix_type,
	  typename
	  std::enable_if<
	    core::meta::is_core_matrix_wrapper<matrix_type>::value &&
	    core::details::traits<matrix_type>::isEigen
	    >::type * = nullptr
	 >
class samplingOperatorIdentity
  : public samplingOperatorBase<samplingOperatorIdentity<matrix_type>>
{
  // this class is for demontration purposes. 
  // later on we need to do sometbing else because
  // we should multiplty by identiy by it makes no sense for perfomance.
  //maybe detect where identity is used and sfinae out

public:
  samplingOperatorIdentity() = default;
  ~samplingOperatorIdentity() = default;
  
private:
  template <typename vector_type>
  void applyImpl(const vector_type & yin,
		 vector_type & yout) const{
    yout = yin;
  }

private:
  friend samplingOperatorBase<samplingOperatorIdentity<matrix_type>>;

};//end class

}//end namespace exp
}//end namespace rom
#endif
