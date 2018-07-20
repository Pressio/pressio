
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_NEWTON_RAPHSON_HPP_
#define SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_NEWTON_RAPHSON_HPP_

#include <memory>

#include "meta/core_meta_static_checks.hpp"

#include "solvers_linear_iterative_factory.hpp"
#include "solvers_nonlinear_iterative_base.hpp"

#include "system_traits.hpp"


namespace solvers {

// Forward declaration
class NonlinearIterativeSolvers;


/**
 * @brief Class that implements a non linear Newton-Raphson solver
 */ 
template<typename SystemT>  
class NonlinearNewtonRaphsonSolver 
  : public NonlinearIterativeSolverBase<
      NonlinearNewtonRaphsonSolver<SystemT>
    >
{

  private:
          
    friend NonlinearIterativeSolvers;
    typedef NonlinearIterativeSolverBase<NonlinearNewtonRaphsonSolver<SystemT>> base_type;

    
  public: 


    NonlinearNewtonRaphsonSolver(NonlinearNewtonRaphsonSolver&& other) : system_(std::move(other.system_)) {}


    template <typename NewSystemT>
    void reassignSystem(const NewSystemT& newSystem) {
      static_assert(details::system_traits<NewSystemT>::is_system, "Error: the object passed to reassignSystem is not a valid system");
      static_assert(details::are_system_compatible<SystemT, NewSystemT>::value, "Error: the new system is not compatible with the existing nonlinear solver");
      system_ = newSystem; 
    }


    template <typename SolverT,
      typename NormT,
      typename VectorT,
      typename std::enable_if<
        core::meta::are_vector_matrix_compatible<VectorT, typename SystemT::matrix_type>::value,
        VectorT
      >::type* = nullptr
    >
    auto solve(const VectorT& xInit) {

      auto dy = system_.residual(xInit);
      auto Ja = system_.jacobian(xInit);

      auto linearSolver = LinearIterativeSolvers::createSolver<SolverT>(Ja);
      linearSolver.setTolerance(this->getLinearSolverTolerance());
      linearSolver.setMaxIterations(this->getLinearSolverMaxIterations());

      int iStep = 1;
      auto xOld = xInit;
      auto xNew = xInit - linearSolver.solve(dy);
 
      while (iStep++ < this->getMaxIterations() && NormT::template compute_norm_difference(xOld, xNew) > this->getTolerance()) {
        xOld = xNew;
        dy = system_.residual(xNew);
        Ja = system_.jacobian(xNew);

        linearSolver.resetLinearSystem(Ja);
        xNew = xNew - linearSolver.solve(dy);
      }

      return xNew;
    }
   
 
    template <typename SolverT,
      typename PrecT,
      typename NormT,
      typename VectorT,
      typename std::enable_if<
        core::meta::are_vector_matrix_compatible<VectorT, typename SystemT::matrix_type>::value,
        VectorT
      >::type* = nullptr
    >
    auto solve(const VectorT& xInit) {

      auto dy = system_.residual(xInit);
      auto Ja = system_.jacobian(xInit);

      auto linearSolver = LinearIterativeSolvers::createSolver<SolverT, PrecT>(Ja);
      linearSolver.setTolerance(this->getLinearSolverTolerance());
      linearSolver.setMaxIterations(this->getLinearSolverMaxIterations());

      int iStep = 1;
      auto xOld = xInit;
      auto xNew = xInit - linearSolver.solve(dy);

      while (iStep++ < this->getMaxIterations() && NormT::template compute_norm_difference(xOld, xNew) > this->getTolerance()) {
        xOld = xNew;
        dy = system_.residual(xNew);
        Ja = system_.jacobian(xNew);

        linearSolver.resetLinearSystem(Ja);
        xNew = xNew - linearSolver.solve(dy);
      }

      return xNew;
    }


  private:

    NonlinearNewtonRaphsonSolver() = delete;
    NonlinearNewtonRaphsonSolver(const SystemT& system) : base_type(), system_(system) {}


  private:
 
    SystemT system_; 
};

} //end namespace solvers

#endif
