
#ifndef SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_
#define SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_

#include <memory>

#include "matrix/core_vector_traits_exp.hpp"
#include "vector/core_vector_traits_exp.hpp"
#include "solvers_linear_iterative_factory.hpp"


namespace solvers {

// Forward declarations
struct NonlinearIterativeSolvers;


/**
 * @brief Class that implements a non linear Newton-Raphson solver
 */ 
template<typename SystemT, 
  typename NormT
>  
class NonlinearNewtonRaphson 
  : public NonlinearIterativeSolverBase<
      NonlinearNewtonRaphson<
        SystemT,
        NormT
      >
    >
{

  private:

    friend NonlinearIterativeSolvers;
    typedef NonlinearIterativeSolverBase<NonlinearNewtonRaphson<SystemT,NormT>> base_type;

    
  public: 

    template <typename NewSystemT>
    void reassignSystem(const NewSystemT& newSystem) {
      static_assert(details::system_traits<T>::is_system, "Error: the object passed to reassignSystem is not a valid system");
      static_assert(meta::are_system_compatible<SystemT, NewSystemT>::value, "Error: the new system is not compatible with the nonlinear solver");
      system_ = newSystem; 
    }


    template <typename SolverT,
      typename VectorT,
      typename std::enable_if<
        typename SystemT::meta::are_vector_matrix_compatible<VectorT, typename SystemT::matrix_type>::value,
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
      auto xNew = linearSolver.solve(dy);
      auto xOld = xInit;

      while (iStep++ < mStep && norm_.check_convergence(xOld, xNew)) {
        xOld = xNew;
        dy = system_.residual(xNew);
        Ja = system_.jacobian(xNew);

        linearSolver.resetLinearSystem(Ja);
        xNew = linearSolver.solve(dy);
      }

      return xNew;
    }
   
 
    template <typename SolverT,
      typename PrecT,
      typename VectorT,
      typename std::enable_if<
        meta::are_vector_matrix_compatible<VectorT, typename SystemT::matrix_type>::value,
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
      auto xNew = linearSolver.solve(dy);
      auto xOld = xInit;

      while (iStep++ < mStep && norm_.check_convergence(xOld, xNew)) {
        xOld = xNew;
        dy = system_.residual(xNew);
        Ja = system_.jacobian(xNew);

        linearSolver.resetLinearSystem(Ja);
        xNew = linearSolver.solve(dy);
      }

      return xNew;
    }


  private:

    NonlinearNewtonRaphson() = delete;
    NonlinearNewtonRaphson(const SystemT& system) : base_type(), system_(system), norm_() {}


  private:
 
    SystemT system_;
    NormT norm_;
};

} //end namespace solvers

#endif
