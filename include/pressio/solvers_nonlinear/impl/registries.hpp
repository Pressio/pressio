
#ifndef SOLVERS_REGISTRIES_HPP_
#define SOLVERS_REGISTRIES_HPP_

#include "levmar_damping.hpp"

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

#define GETMETHOD(N) \
  template<class Tag, std::enable_if_t< std::is_same<Tag, Tag##N >::value, int> = 0> \
  auto & get(){ return d##N##_; } \
  template<class Tag, std::enable_if_t< std::is_same<Tag, Tag##N >::value, int> = 0> \
  const auto & get() const { return d##N##_; }


template<class SystemType, class InnSolverType>
class RegistryNewton
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;

  using Tag1 = nonlinearsolvers::CorrectionTag;
  using Tag2 = nonlinearsolvers::InitialGuessTag;
  using Tag3 = nonlinearsolvers::ResidualTag;
  using Tag4 = nonlinearsolvers::JacobianTag;
  using Tag5 = nonlinearsolvers::InnerSolverTag;
  using Tag6 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  state_t d2_;
  r_t d3_;
  j_t d4_;
  utils::InstanceOrReferenceWrapper<InnSolverType> d5_;
  SystemType const * d6_;

public:
  template<class _InnSolverType>
  RegistryNewton(const SystemType & system, _InnSolverType && innS)
    : d1_(system.createState()),
      d2_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(std::forward<_InnSolverType>(innS)),
      d6_(&system){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	   Tag1, Tag2, Tag3, Tag4, Tag5, Tag6>::value) < 6;
  }

  GETMETHOD(1)
  GETMETHOD(2)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
};

template<class SystemType, class InnSolverType>
class RegistryGaussNewtonNormalEqs
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;

  using Tag1 = nonlinearsolvers::CorrectionTag;
  using Tag2 = nonlinearsolvers::InitialGuessTag;
  using Tag3 = nonlinearsolvers::ResidualTag;
  using Tag4 = nonlinearsolvers::JacobianTag;
  using Tag5 = nonlinearsolvers::GradientTag;
  using Tag6 = nonlinearsolvers::HessianTag;
  using Tag7 = nonlinearsolvers::InnerSolverTag;
  using Tag8 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  state_t d2_;
  r_t d3_;
  j_t d4_;
  gradient_t d5_;
  hessian_t d6_;
  utils::InstanceOrReferenceWrapper<InnSolverType> d7_;
  SystemType const * d8_;

public:
  template<class _InnSolverType>
  RegistryGaussNewtonNormalEqs(const SystemType & system, _InnSolverType && innS)
    : d1_(system.createState()),
      d2_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(system.createState()),
      d6_( hg_default::createHessian(system.createState()) ),
      d7_(std::forward<_InnSolverType>(innS)),
      d8_(&system){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	   Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8>::value) < 8;
  }

  GETMETHOD(1)
  GETMETHOD(2)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
};

template<class SystemType, class InnSolverType, class WeightingOpType>
class RegistryWeightedGaussNewtonNormalEqs
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;

  using Tag1  = nonlinearsolvers::CorrectionTag;
  using Tag2  = nonlinearsolvers::InitialGuessTag;
  using Tag3  = nonlinearsolvers::ResidualTag;
  using Tag4  = nonlinearsolvers::JacobianTag;
  using Tag5  = nonlinearsolvers::WeightedResidualTag;
  using Tag6  = nonlinearsolvers::WeightedJacobianTag;
  using Tag7  = nonlinearsolvers::GradientTag;
  using Tag8  = nonlinearsolvers::HessianTag;
  using Tag9  = nonlinearsolvers::InnerSolverTag;
  using Tag10 = nonlinearsolvers::WeightingOperatorTag;
  using Tag11 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  state_t d2_;
  r_t d3_;
  j_t d4_;
  r_t d5_;
  j_t d6_;
  gradient_t d7_;
  hessian_t d8_;
  utils::InstanceOrReferenceWrapper<InnSolverType> d9_;
  utils::InstanceOrReferenceWrapper<WeightingOpType> d10_;
  SystemType const * d11_;

public:
  template<class _InnSolverType, class _WeightingOpType>
  RegistryWeightedGaussNewtonNormalEqs(const SystemType & system,
				       _InnSolverType && innS,
				       _WeightingOpType && weigher)
    : d1_(system.createState()),
      d2_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(system.createResidual()),
      d6_(system.createJacobian()),
      d7_(system.createState()),
      d8_( hg_default::createHessian(system.createState()) ),
      d9_(std::forward<InnSolverType>(innS)),
      d10_(std::forward<_WeightingOpType>(weigher)),
      d11_(&system){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	    Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8, Tag9, Tag10, Tag11>::value) < 11;
  }

  GETMETHOD(1)
  GETMETHOD(2)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
  GETMETHOD(9)
  GETMETHOD(10)
  GETMETHOD(11)
};


template<class SystemType, class QRSolverType>
class RegistryGaussNewtonQr
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using QTr_t      = state_t; // type of Q^T*r
  using gradient_t = state_t; // type of J^T r

  using Tag1 = nonlinearsolvers::CorrectionTag;
  using Tag2 = nonlinearsolvers::InitialGuessTag;
  using Tag3 = nonlinearsolvers::ResidualTag;
  using Tag4 = nonlinearsolvers::JacobianTag;
  using Tag5 = nonlinearsolvers::GradientTag;
  using Tag6 = nonlinearsolvers::impl::QTransposeResidualTag;
  using Tag7 = nonlinearsolvers::InnerSolverTag;
  using Tag8 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  state_t d2_;
  r_t d3_;
  j_t d4_;
  gradient_t d5_;
  QTr_t d6_;
  utils::InstanceOrReferenceWrapper<QRSolverType> d7_;
  SystemType const * d8_;

public:
  template<class QRType>
  RegistryGaussNewtonQr(const SystemType & system, QRType && qrs)
    : d1_(system.createState()),
      d2_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(system.createState()),
      d6_(system.createState()),
      d7_(std::forward<QRType>(qrs)),
      d8_(&system){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	   Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8>::value) < 8;
  }

  GETMETHOD(1)
  GETMETHOD(2)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
};

template<class SystemType, class InnSolverType>
class RegistryLevMarNormalEqs
{
  using scalar_t   = scalar_of_t<SystemType>;
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;
  using lm_damp_t  = LevenbergMarquardtDamping<scalar_t>;

  using Tag1 = nonlinearsolvers::CorrectionTag;
  using Tag2 = nonlinearsolvers::InitialGuessTag;
  using Tag3 = nonlinearsolvers::ResidualTag;
  using Tag4 = nonlinearsolvers::JacobianTag;
  using Tag5 = nonlinearsolvers::GradientTag;
  using Tag6 = nonlinearsolvers::HessianTag;
  using Tag7 = nonlinearsolvers::LevenbergMarquardtUndampedHessianTag;
  using Tag8 = nonlinearsolvers::LevenbergMarquardtDampingTag;
  using Tag9 = nonlinearsolvers::InnerSolverTag;
  using Tag10 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  state_t d2_;
  r_t d3_;
  j_t d4_;
  gradient_t d5_;
  hessian_t d6_;
  hessian_t d7_;
  lm_damp_t d8_;
  utils::InstanceOrReferenceWrapper<InnSolverType> d9_;
  SystemType const * d10_;

public:
  template<class _InnSolverType>
  RegistryLevMarNormalEqs(const SystemType & system, _InnSolverType && innS)
    : d1_(system.createState()),
      d2_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(system.createState()),
      d6_( hg_default::createHessian(system.createState()) ),
      d7_( hg_default::createHessian(system.createState()) ),
      d8_{},
      d9_(std::forward<_InnSolverType>(innS)),
      d10_(&system){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	    Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8, Tag9, Tag10>::value) < 10;
  }

  GETMETHOD(1)
  GETMETHOD(2)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
  GETMETHOD(9)
  GETMETHOD(10)
};

}}}
#endif
