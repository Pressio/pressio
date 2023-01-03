
#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <class HessianType, class GradientType, class ResidualNormValueType>
class HessianGradientOperators{
public:
  using residual_norm_value_type = ResidualNormValueType;

private:
  GradientType g_;
  HessianType H_;
  ResidualNormValueType residualNorm_;

public:
  template <typename SystemType>
  HessianGradientOperators(const SystemType & system)
    : g_(system.createGradient()),
      H_(system.createHessian())
  {
    ::pressio::ops::set_zero(g_);
    ::pressio::ops::set_zero(H_);
  }

public:
  const HessianType & hessianCRef() const { return H_; }
  const GradientType & gradientCRef() const { return g_; }

  template< typename SystemType, typename StateType>
  void residualNorm(const SystemType & system,
		    const StateType & state,
		    ResidualNormValueType & residualNorm) const
  {
    system.residualNorm(state, ::pressio::Norm::L2, residualNorm);
  }

  template<typename SystemType, typename StateType>
  void compute(const SystemType & sys,
	       const StateType & state,
	       bool recomputeSystemJacobian = true)
  {
    if (recomputeSystemJacobian){
      sys.hessian(state, H_);
    }

    sys.gradient(state, g_, ::pressio::Norm::L2,
		 residualNorm_, recomputeSystemJacobian);

    // scale because of sign convention
    using g_scalar_type = typename ::pressio::Traits<GradientType>::scalar_type;
    ::pressio::ops::scale(g_, ::pressio::utils::Constants<g_scalar_type>::negOne());
  }
};

template <class ResidualType, class JacobianType>
class ResidualJacobianOperators
{
public:
  using residual_norm_value_type = decltype(::pressio::ops::norm2(std::declval<ResidualType>()));
  using jacobian_type = JacobianType;
  using residual_type = ResidualType;

private:
  ResidualType r_;
  JacobianType J_;
  // auxR_ is used for residualNorm method so that we don't modify
  // the real operator r_ which must be the same once computeOperators is called
  mutable ResidualType auxR_;
  residual_norm_value_type residualNorm_;

public:
  template <class SystemType>
  ResidualJacobianOperators(const SystemType & system)
    : r_(system.createResidual()),
      J_(system.createJacobian()),
      auxR_(::pressio::ops::clone(r_))
  {
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(J_);
  }

public:
  const ResidualType & residualCRef() const { return r_; }
  const JacobianType & jacobianCRef() const { return J_; }

  template<class SystemType, class StateType>
  void compute(const SystemType & sys,
	       const StateType & state,
	       bool recomputeSystemJacobian=true)
  {
    sys.residual(state, r_);
    residualNorm_ = ::pressio::ops::norm2(r_);
    if (std::isnan(residualNorm_)){
      throw ::pressio::eh::ResidualHasNans();
    }

    if (recomputeSystemJacobian){
      sys.jacobian(state, J_);
    }
  }

  residual_norm_value_type residualNorm() const{
    return residualNorm_;
  }

  template< class SystemType, class state_t>
  void residualNorm(const SystemType & system,
		    const state_t & state,
		    residual_norm_value_type & residualNorm) const
  {
    system.residual(state, auxR_);
    residualNorm = ::pressio::ops::norm2(auxR_);

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_HPP_
