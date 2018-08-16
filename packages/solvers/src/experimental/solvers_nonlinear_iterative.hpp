
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_HPP

#include <type_traits>

#include "solvers_nonlinear_base.hpp"
#include "solvers_nonlinear_factory.hpp"


namespace solvers {


struct NonLinearSolvers; // Fwd declaration


template <
  typename PolicyT
>
class NonLinearIterativeSolver 
  : public NonLinearSolverBase<
      NonLinearIterativeSolver<
        PolicyT
      >
    >
{

  private:

    friend NonLinearSolvers;
  	typedef NonLinearSolverBase<NonLinearIterativeSolver<PolicyT>> base_type;


  public:


    /**
     *
     *
     */
    template <
      typename SolverT,
      typename PrecT,
      typename NormT,
      typename SystemT,
      typename VectorT
    >
    auto solve_(const SystemT& sys, const VectorT& b) {
      return 0;
    }


    /**
     * @brief Get the maximum number of iterations of the non linear iterative solver
     *
     */
    uint getMaxNonLinearIterations() {
      return maxNonLinearIterations_;
    }
  

    /**
     * @brief Get the tolerance of the non linear iterative solver
     *
     */
    double getNonLinearTolerance() {
      return nonLinearTolerance_;
    }
    

    /**
     * @brief Set the maximum number of iterations of the non linear iterative solver
     *
     * @param maxNonLinearIterations represents number of iterations of the non linear solver
     */
    void setMaxNonLinearIterations(uint maxNonLinearIterations) {
      maxNonLinearIterations_ = maxNonLinearIterations;
    }

    /**
     * @brief Set the tolerance of the non linear iterative solver
     *
     * @param tolerance is the tolerance of the non linear solver
     */
    void setNonLinearTolerance(double nonLinearTolerance) {
      nonLinearTolerance_ = abs(nonLinearTolerance);
    }
  	

  protected:

  	NonLinearIterativeSolver() : base_type(), maxNonLinearIterations_(100), nonLinearTolerance_(1.0e-5) {}


  private:

    uint maxNonLinearIterations_;
    double nonLinearTolerance_;
};
	
} // end namespace solvers

#endif