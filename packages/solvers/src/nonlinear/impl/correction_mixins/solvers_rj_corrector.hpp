
#ifndef PRESSIO_SOLVERS_RJ_CORRECTOR_IMPL_HPP_
#define PRESSIO_SOLVERS_RJ_CORRECTOR_IMPL_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<
  typename T, typename state_t, typename lin_solver_t, ::pressio::Norm normType
  >
class RJCorrector : public T
{
  using sc_t = typename ::pressio::containers::details::traits<state_t>::scalar_t;

  state_t correction_ = {};
  lin_solver_t & solverObj_;
  sc_t residualNorm_ = {};

public:
  static constexpr auto normType_ = normType;

  RJCorrector() = delete;

  template <typename system_t>
  RJCorrector(const system_t & system, 
        const state_t & state, 
        lin_solver_t & solverObj)
    : T(system, state), correction_(state), solverObj_(solverObj){}

public:
  template <typename system_t>
  void computeCorrection(const system_t & sys, state_t & state)
  {
    T::computeOperators(sys, state, normType, residualNorm_);

    auto & r = T::getResidual();
    auto & J = T::getJacobian();

    // solve: R correction = Q^T Residual
    solverObj_.solve(J, r, correction_);
    // scale by -1 for sign convention
    pressio::ops::scale(correction_, utils::constants<sc_t>::negOne() );

    std::cout << std::fixed
    	      << std::setprecision(15)
    	      << residualNorm_ << " "
    	      << pressio::ops::norm2(correction_)
    	      << std::endl;
  }

  const state_t & viewCorrection() const{ return correction_; }
  const sc_t residualNormCurrentCorrectionStep() const{ return residualNorm_; }

  template< typename system_t>
  void residualNorm(const system_t & system, const state_t & state, sc_t & result){
    T::residualNorm(system, state, normType, result);
  }

};

}}}}
#endif
