
#ifndef ODE_IMPLICIT_RESIDUAL_POLICY_BASE_HPP_
#define ODE_IMPLICIT_RESIDUAL_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{   

//-----------------------------------------------------------------

template <int numAuxStates, int numAuxRHS, typename derived_t,
	  typename enable = void>
class implicitResidualPolicyBase;

//-----------------------------------------------------------------
  
template <int numAuxStates, int numAuxRHS, typename derived_t>
class implicitResidualPolicyBase<numAuxStates, numAuxRHS, derived_t,
				 typename std::enable_if<
				   numAuxRHS==0 >::type>
{
public:
  template <typename state_type,
	    typename residual_type,
	    typename model_type,
	    typename time_type>
  void compute(const state_type & y,
	      residual_type & R,
	      const std::array<state_type, numAuxStates> & auxYs,
	      model_type & model,
	      time_type t,
	      time_type dt)
  {
    this->underlying().computeImpl(y, R, auxYs, model, t, dt);
  }

private:
  friend derived_t;

  implicitResidualPolicyBase() = default;
  ~implicitResidualPolicyBase() = default;
  
  derived_t & underlying(){
    return static_cast<derived_t& >( *this );
  }
  derived_t const & underlying() const{
    return static_cast<derived_t const & >( *this );
  } 
};//end class

//-----------------------------------------------------------------


template <int numAuxStates, int numAuxRHS, typename derived_t>
class implicitResidualPolicyBase<numAuxStates, numAuxRHS, derived_t,
				 typename std::enable_if<
				   numAuxStates!=0 &&
				   numAuxRHS!=0 >::type>
{
public:  
  template <typename state_type, typename residual_type,
  	    typename model_type, typename time_type>
  void compute(const state_type & y,
  	       residual_type & R,
  	       const std::array<state_type, numAuxStates> & auxYs,
  	       const std::array<residual_type, numAuxRHS> & auxRHSs,
  	       model_type & model,
  	       time_type t,
  	       time_type dt)
  {
    this->underlying().computeImpl(y, R, auxYs, auxRHSs, model, t, dt);
  }
  
private:
  friend derived_t;

  implicitResidualPolicyBase() = default;
  ~implicitResidualPolicyBase() = default;
  
  derived_t & underlying(){
    return static_cast<derived_t& >( *this );
  }
  derived_t const & underlying() const{
    return static_cast<derived_t const & >( *this );
  } 
};//end class
//-----------------------------------------------------------------







  
}//end namespace polices
}//end namespace ode  
#endif
