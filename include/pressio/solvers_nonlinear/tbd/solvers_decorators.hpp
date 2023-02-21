
#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <class HessianType, class GradientType, class ParentType>
class HessianGradientDecorator : public ParentType
{
private:
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();

  GradientType g_;
  HessianType H_;

public:
  template <typename SystemType>
  HessianGradientDecorator(const SystemType & system)
    : ParentType(system),
      g_(system.createState()),
      H_(::pressio::ops::product<HessianType>
	 (pT, pnT, ::pressio::utils::Constants<
	  typename ::pressio::Traits<typename ParentType::jacobian_type>::scalar_type>::one(),
	  ParentType::jacobianCRef()))
  {
    ::pressio::ops::set_zero(g_);
    ::pressio::ops::set_zero(H_);
  }

public:
  const HessianType & hessianCRef() const   { return H_; }
  const GradientType & gradientCRef() const { return g_; }

public:
  template<class SystemType, class StateType>
  void compute(const SystemType & systemObj,
	       const StateType & state,
	       bool recomputeSystemJacobian = true)
  {
    ParentType::compute(systemObj, state, recomputeSystemJacobian);
    // gradient always computed because residual always changes
    this->_computeGradient();

    if (recomputeSystemJacobian){
      this->_computeHessian();
    }
  }

private:
  void _computeHessian()
  {
    using H_scalar_type = typename ::pressio::Traits<HessianType>::scalar_type;
    constexpr auto beta  = ::pressio::utils::Constants<H_scalar_type>::zero();
    using j_t = typename ParentType::jacobian_type;
    using J_scalar_type = typename ::pressio::Traits<j_t>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<J_scalar_type>::one();
    ::pressio::ops::product(pT, pnT, alpha, ParentType::jacobianCRef(), beta, H_);
  }

  void _computeGradient()
  {
    using g_scalar_type = typename ::pressio::Traits<GradientType>::scalar_type;
    constexpr auto beta  = ::pressio::utils::Constants<g_scalar_type>::zero();

    using j_t = typename ParentType::jacobian_type;
    using J_scalar_type = typename ::pressio::Traits<j_t>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<J_scalar_type>::one();
    // compute gradient (g_ = J^T r)
    ::pressio::ops::product(pT, alpha, ParentType::jacobianCRef(), ParentType::residualCRef(), beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::Constants<g_scalar_type>::negOne());
  }
};


template <class HessianType, class GradientType, class WeightFuncType, class ParentType>
class WeightedHessianGradientDecorator : public ParentType
{
  using residual_type = typename ParentType::ResidualType;
  using jacobian_type = typename ParentType::JacobianType;

private:
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();

  GradientType g_;
  HessianType H_;
  ::pressio::utils::InstanceOrReferenceWrapper<WeightFuncType> functorW_;
  residual_type Wr_;
  jacobian_type WJ_;

public:
  template <typename SystemType, class WType>
  WeightedHessianGradientDecorator(const SystemType & system,
				   WType && w)
    : ParentType(system),
      g_(system.createState()),
      H_(::pressio::ops::product<HessianType>
	 (pT, pnT, ::pressio::utils::Constants<
	  typename ::pressio::Traits<typename ParentType::jacobian_type>::scalar_type>::one(),
	  ParentType::jacobianCRef())),
      functorW_(std::forward<WType>(w)),
      Wr_(system.createResidual()),
      WJ_(system.createJacobian())
  {
    ::pressio::ops::set_zero(g_);
    ::pressio::ops::set_zero(H_);
    ::pressio::ops::set_zero(Wr_);
    ::pressio::ops::set_zero(WJ_);
  }

public:
  const HessianType & hessianCRef() const   { return H_; }
  const GradientType & gradientCRef() const { return g_; }

public:
  template<class SystemType, class StateType>
  void compute(const SystemType & systemObj,
	       const StateType & state,
	       bool recomputeSystemJacobian = true)
  {
    ParentType::compute(systemObj, state, recomputeSystemJacobian);
    functorW_(ParentType::residualCRef(), Wr_);

    // gradient always computed because residual always changes
    this->_computeGradient();

    if (recomputeSystemJacobian){
      functorW_(ParentType::jacobianCRef(), WJ_);
      this->_computeHessian();
    }
  }

private:
  void _computeHessian()
  {
    using H_scalar_type = typename ::pressio::Traits<HessianType>::scalar_type;
    constexpr auto beta  = ::pressio::utils::Constants<H_scalar_type>::zero();

    using J_scalar_type = typename ::pressio::Traits<jacobian_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<J_scalar_type>::one();
    ::pressio::ops::product(pT, pnT, alpha,
			    ParentType::jacobianCRef(), WJ_,
			    beta, H_);
  }

  void _computeGradient()
  {
    using g_scalar_type = typename ::pressio::Traits<GradientType>::scalar_type;
    constexpr auto beta  = ::pressio::utils::Constants<g_scalar_type>::zero();

    using J_scalar_type = typename ::pressio::Traits<jacobian_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<J_scalar_type>::one();
    // compute gradient (g_ = J^T W r)
    ::pressio::ops::product(pT, alpha,
			    ParentType::jacobianCRef(), Wr_, beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::Constants<g_scalar_type>::negOne());
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_
