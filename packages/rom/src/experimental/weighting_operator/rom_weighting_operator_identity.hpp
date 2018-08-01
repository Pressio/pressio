
#ifndef ROM_WEIGHTING_OPERATOR_IDENTITY_HPP_
#define ROM_WEIGHTING_OPERATOR_IDENTITY_HPP_

#include "rom_ConfigDefs.hpp"
#include "./rom_weighting_operator_base.hpp"
#include <type_traits>

namespace rom{
namespace experimental{

template<typename matrix_type,
	 typename phi_op_type,
	  typename
	  std::enable_if<
	    core::meta::is_core_matrix_wrapper<matrix_type>::value &&
	    core::details::traits<matrix_type>::isEigen
	    >::type * = nullptr
	 >
class weightingOperatorIdentity
  : public weightingOperatorBase<
  weightingOperatorIdentity<matrix_type, phi_op_type>>
{
public:
  weightingOperatorIdentity(phi_op_type & phiObj,
			    size_t nr, size_t nc)
    : phiOp_(&phiObj), tmpM_(nr, nc)
  {}

  ~weightingOperatorIdentity() = default;
  
private:
  template <typename T>
  void leftMultiplyImpl(const T & yin, T & yout) 
  {
    phiOp_->project(yin, yout);
  }

  template <typename T>
  void leftRightMultiplyImpl(const T & Min, T & Mout) 
  {
    // static_assert(core::details::traits<matrix_type>::isSparse &&
    // 		  core::details::traits<matrix_type>::isEigen, "");
    phiOp_->rightMultiply(Min, tmpM_);
    phiOp_->project(tmpM_, Mout);
  }
  
private:
  friend weightingOperatorBase<
  weightingOperatorIdentity<matrix_type, phi_op_type> >;

private:
  phi_op_type * phiOp_;
  matrix_type tmpM_;
  
};//end class

}//end namespace exp
}//end namespace rom
#endif
