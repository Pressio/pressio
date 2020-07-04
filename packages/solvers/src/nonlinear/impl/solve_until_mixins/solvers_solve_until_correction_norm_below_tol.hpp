

#ifndef PRESSIO_SOLVERS_STOP_WHEN_CORRECTION_NORM_BELOW_TOL_HPP_
#define PRESSIO_SOLVERS_STOP_WHEN_CORRECTION_NORM_BELOW_TOL_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<typename sc_t, typename T>
class SolveUntilCorrectionNormBelowTol
  : public T,
    public IterativeBase< SolveUntilCorrectionNormBelowTol<sc_t, T>, sc_t >
{
  using this_t = SolveUntilCorrectionNormBelowTol<sc_t, T>;
  using iterative_base_t = IterativeBase<this_t, sc_t>;
  using typename iterative_base_t::iteration_t;

  iteration_t iStep_ = {};

public:
  SolveUntilCorrectionNormBelowTol() = delete;

  template <typename ...Args>
  SolveUntilCorrectionNormBelowTol(Args &&... args)
    : T(std::forward<Args>(args)...){}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    iStep_ = 0;
    while (++iStep_ <= iterative_base_t::maxIters_)
    {
      T::computeCorrection(sys, state);
      T::updateState(sys, state);

      const auto & correction = T::viewCorrection();
      const auto correctionNorm = ::pressio::ops::norm2(correction);
      if (correctionNorm < iterative_base_t::tolerance_)
        break;
    }
  }

private:
  iteration_t getNumIterationsExecutedImpl() const {
    return iStep_;
  }  

};

}}}}
#endif
