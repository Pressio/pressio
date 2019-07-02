
#ifndef ROM_LSPG_STEADY_SYSTEM_HPP_
#define ROM_LSPG_STEADY_SYSTEM_HPP_

namespace pressio{ namespace rom{

template<
  typename app_type,
  typename lspg_state_type,
  typename lspg_residual_type,
  typename lspg_jacobian_type,
  typename residual_policy_type,
  typename jacobian_policy_type
  >
class LSPGSteadySystem<
  app_type, lspg_state_type, lspg_residual_type,
  lspg_jacobian_type, residual_policy_type,
  jacobian_policy_type
  >{

  const app_type & app_;
  const residual_policy_type & residualEvaluator_;
  const jacobian_policy_type & jacobianEvaluator_;

public:
  // these need to be public because are detected by solver
  using scalar_type = typename app_type::scalar_type;
  using state_type	= lspg_state_type;
  using residual_type	= lspg_residual_type;
  using jacobian_type	= lspg_jacobian_type;

public:
  LSPGSteadySystem() = delete;
  ~LSPGSteadySystem() = default;

  LSPGSteadySystem(const app_type & appIn,
		   const residual_policy_type & resPolicyObj,
		   const jacobian_policy_type & jacPolicyObj)
    : app_(appIn),
      residualEvaluator_(resPolicyObj),
      jacobianEvaluator_(jacPolicyObj){}

public:
  void residual(const lspg_state_type & y,
		lspg_residual_type & R) const{
    (this->residualEvaluator_).template operator()(y, R, app_);
  }

  void jacobian(const lspg_state_type & y,
		lspg_jacobian_type & J) const{
    (this->jacobianEvaluator_).template operator()(y, J, app_);
  }

  lspg_residual_type residual(const lspg_state_type & y) const{
    return (this->residualEvaluator_).template operator()(y, app_);
  }

  lspg_jacobian_type jacobian(const lspg_state_type & y) const{
    return (this->jacobianEvaluator_).template operator()(y, app_);
  }
};//end class

}} // end namespace pressio::ode
#endif
