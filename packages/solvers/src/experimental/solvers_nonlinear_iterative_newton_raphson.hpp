
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_NEWTON_RAPHSON_HPP_
#define SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_NEWTON_RAPHSON_HPP_

#include <memory>

#include "meta/core_meta_static_checks.hpp"

#include "solvers_linear_iterative_factory.hpp"
#include "solvers_nonlinear_iterative_base.hpp"

#include "solvers_traits.hpp"
#include "system_traits.hpp"


namespace solvers {

// Forward declaration
class NonlinearIterativeSolvers;


/**
 * @brief Class that implements a non linear Newton-Raphson solver
 */  
class NonlinearNewtonRaphsonSolver 
  : public NonlinearIterativeSolverBase<NonlinearNewtonRaphsonSolver>
{

  private:
          
    friend NonlinearIterativeSolvers;
    typedef NonlinearIterativeSolverBase<NonlinearNewtonRaphsonSolver> base_type;

    
  public: 

    /// Move ctor
    NonlinearNewtonRaphsonSolver(NonlinearNewtonRaphsonSolver&& other) {};
    

    /// Uses default preconditioner and L2 norm 
    template <typename SolverT,
      typename SystemT,
      typename VectorT,
      typename std::enable_if<
        details::system_traits<SystemT>::is_system &&
        core::meta::are_vector_matrix_compatible<VectorT, typename SystemT::matrix_type>::value,
        VectorT
      >::type* = nullptr
    >
    auto solve(const SystemT& system, const VectorT& xInit) {
      return this->template solve<SolverT, typename linear::DefaultPreconditioner, L2Norm>(system, xInit);
    }

       
    // Specifies preconditioner and norm to be used
    template <typename SolverT,
      typename PrecT,
      typename NormT,
      typename SystemT,
      typename VectorT,
      typename std::enable_if<
      details::system_traits<SystemT>::is_system &&
        core::meta::are_vector_matrix_compatible<VectorT, typename SystemT::matrix_type>::value,
        VectorT
      >::type* = nullptr
    >
    auto solve(const SystemT& system, const VectorT& xInit) {

      auto dy = system.residual(xInit);
      auto Ja = system.jacobian(xInit);

      auto linearSolver = LinearIterativeSolvers::createSolver<SolverT, PrecT>(Ja);
      linearSolver.setTolerance(this->getLinearSolverTolerance());
      linearSolver.setMaxIterations(this->getLinearSolverMaxIterations());

      int iStep = 1;
      auto xOld = xInit;
      auto xNew = xInit - linearSolver.solve(dy);

      while (iStep++ < this->getMaxIterations() && NormT::template compute_norm_difference(xOld, xNew) > this->getTolerance()) {
        xOld = xNew;
        dy = system.residual(xNew);
        Ja = system.jacobian(xNew);

        linearSolver.resetLinearSystem(Ja);
        xNew = xNew - linearSolver.solve(dy);
      }

      return xNew;
    }


  private:

    NonlinearNewtonRaphsonSolver() : base_type() {};

};

} //end namespace solvers

#endif
