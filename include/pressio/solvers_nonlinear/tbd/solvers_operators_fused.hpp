
#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_FUSED_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_FUSED_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <class HessianType, class GradientType, class ResidualNormValueType>
class FusedHessianGradientOperators{
public:
  using residual_norm_value_type = ResidualNormValueType;

private:
  GradientType g_;
  HessianType H_;
  ResidualNormValueType residualNorm_;

public:
  template <typename SystemType>
  FusedHessianGradientOperators(const SystemType & system)
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
    sys.hessianAndGradient(state, H_, g_,
			   ::pressio::Norm::L2,
			   residualNorm_, recomputeSystemJacobian);

    // scale because of sign convention
    using g_scalar_type = typename ::pressio::Traits<GradientType>::scalar_type;
    ::pressio::ops::scale(g_, ::pressio::utils::Constants<g_scalar_type>::negOne());
  }
};

template <class ResidualType, class JacobianType>
class FusedResidualJacobianOperators
{
public:
  using residual_norm_value_type = decltype(::pressio::ops::norm2(std::declval<ResidualType>()));

private:
  // J_ is mutable because for the fused_res_jac api, to compute the
  // residualNorm (which is a const method), we call
  // residualAndJacobian without recomputing J_ but J_ still needs to be passed
  ResidualType r_;
  mutable JacobianType J_;
  mutable ResidualType auxR_;
  residual_norm_value_type residualNorm_;

public:
  template <class SystemType>
  FusedResidualJacobianOperators(const SystemType & system)
    : r_( system.createResidual() ),
      J_( system.createJacobian() ),
      auxR_(::pressio::ops::clone(r_))
  {
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(J_);
  }

public:
  const ResidualType & residualCRef() const { return r_; }
  const JacobianType & jacobianCRef() const { return J_; }

  template<class SystemType, class state_t>
  void compute(const SystemType & sys,
	       const state_t & state,
	       bool recomputeSystemJacobian=true)
  {
    sys.residualAndJacobian(state, r_, J_, recomputeSystemJacobian);
    residualNorm_ = ::pressio::ops::norm2(r_);

    if (std::isnan(residualNorm_)){
      throw ::pressio::eh::ResidualHasNans();
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
    system.residualAndJacobian(state, auxR_, J_, false);
    residualNorm = ::pressio::ops::norm2(auxR_);

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_HPP_
