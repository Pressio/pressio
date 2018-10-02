
#ifndef CORE_SHARED_BASE_IN_PLACE_OPS_HPP_
#define CORE_SHARED_BASE_IN_PLACE_OPS_HPP_

#include "../core_vector_traits.hpp"
  
namespace rompp{
namespace core{
    
template<typename derived_type>
class ContainerInPlaceOpsBase
  : private core::details::CrtpBase<
  ContainerInPlaceOpsBase<derived_type>>{

  using sc_t = typename details::traits<derived_type>::scalar_t;

public:
  template <typename op_t,
	    typename T = derived_type>
  void inPlaceOp(sc_t a1, sc_t a2,
		 const T & other){
    // this = a1*this op a2*other;
    this->underlying().template inPlaceOpImpl<op_t,T>(a1, a2, other);
  }

  template <typename op_t,
  	    typename T = derived_type>
  void inPlaceOp(sc_t a1, const T & x1,
  		 sc_t a2, const T & x2){
    // this = a1*x1 op a2*x2;
    this->underlying().template inPlaceOpImpl<op_t,T>(a1, x1, a2, x2);
  }
  
  template <typename op0_t,
	    typename T = derived_type,
	    typename op1_t = op0_t,
  	    typename op2_t = op0_t,
  	    typename op3_t = op0_t>
  void inPlaceOp(sc_t a0,
  		 sc_t a1, const T & x1,
  		 sc_t a2, const T & x2,
  		 sc_t a3, const T & x3,
  		 sc_t a4, const T & x4){
    // this = a0 * this op0 (a1*x1) op1 (a2*x2) op2 (a3*x3) op3 (a4*x4)
    this->underlying().template inPlaceOpImpl<op0_t, T, op1_t,
  					      op2_t, op3_t>(a0, a1, x1,
							    a2, x2,
							    a3, x3,
							    a4, x4);
  }
  
private:
  friend derived_type;
  friend core::details::CrtpBase<ContainerInPlaceOpsBase<derived_type>>;

  VectorMathBase() = default;
  ~VectorMathBase() = default;

};//end class

  
} // end namespace core
}//end namespace rompp
#endif
