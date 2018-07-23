
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_BASE_HPP_
#define SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_BASE_HPP_


namespace solvers {


/**
 * @brief Base class for nonlinear iterative solver using CRTP
 *
 * @section DESCRIPTION
 *
 * This class defines the methods common to different nonlinear iterative solvers.
 *  
 * Objects of the class cannot be created directly. To create a solver,
 * use the factory class NonlinearIterativeSolvers.
 */
template <typename Derived>
class NonlinearIterativeSolverBase {

  public: 

    // Uses default preconditioner and L2 norm
    template <typename SolverT,
      typename SystemT,
      typename VectorT
    > 
    auto solve(const SystemT& system, const VectorT& xInit) {
      return this->underlying().template solve<SolverT>(system, xInit);
    }


    // Specifies linear solver, preconditioner and norm to be used
    template <typename SolverT,
      typename PrecT,
      typename NormT,
      typename SystemT,
      typename VectorT,
      typename std::enable_if<
        core::meta::are_vector_matrix_compatible<VectorT, typename SystemT::matrix_type>::value,
        VectorT
      >::type* = nullptr
    >
    auto solve(const SystemT& system, const VectorT& xInit) {
      return this->underlying().template solve<SolverT, PrecT, NormT>(system, xInit);
    }


    int getMaxIterations() {
      return maxIters_;
    }


    int getLinearSolverMaxIterations() {
      return linearSolverMaxIters_;
    }


    void setMaxIterations(int maxIters) {
      maxIters_ = maxIters;
    }


    void setLinearSolverMaxIterations(int maxIters) {
      linearSolverMaxIters_ = maxIters;
    }


    double getTolerance() {
      return tolerance_;
    }

    inline double getLinearSolverTolerance() {
      return linearSolverTolerance_;
    }


    void setTolerance(double tolerance) {
      tolerance_ = tolerance;
    }


    void setLinearSolverTolerance(double tolerance) {
      linearSolverTolerance_ = tolerance;
    }


  protected:

    NonlinearIterativeSolverBase() : maxIters_(100), linearSolverMaxIters_(100), tolerance_(1.0e-5), linearSolverTolerance_(1.0e-7) {}
    ~NonlinearIterativeSolverBase() = default;

    NonlinearIterativeSolverBase(const NonlinearIterativeSolverBase&) = delete;


  private:

    Derived& underlying() {
      return static_cast<Derived&>(*this);
    }
  

    Derived const& underlying() const {
      return static_cast<Derived const&>(*this);
    }


  private:

    int maxIters_, linearSolverMaxIters_;
    double tolerance_, linearSolverTolerance_;
      
};

} //end namespace solvers

#endif
