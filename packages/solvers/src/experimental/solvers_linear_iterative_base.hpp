
#ifndef SOLVERS_LINEAR_ITERATIVE_BASE_HPP_
#define SOLVERS_LINEAR_ITERATIVE_BASE_HPP_

namespace solvers {


/**
 * @brief Base class for linear iterative solver implemented through CRTP
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a linear iterative solver. 
 * Objects of the class cannot be created directly. To create a solver,
 * use the factory class LinearIterativeSolvers.
 */
template<typename Derived>
class LinearIterativeSolverBase {

<<<<<<< HEAD
=======
  private:

    typedef Derived::matrix_type matrix_type;


>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
  public: 


    /**
     * @brief  Initialize a new linear solver
     *
     * @param  A Matrix representing the linear system to solve
     */
<<<<<<< HEAD
    template <typename T>
    void resetLinearSystem(const T& A) {
=======
    void resetLinearSystem(const matrix_type& A) {
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
      this->underlying().resetLinearSystem(A);
    }


    /**
     * @brief  Solve the linear system
     *
     * @param  B is the RHS vector
     * @return Solution vector
     */
    template <typename T>
    auto solve(const T& b) {
      return this->underlying().solve(b);  
    }


    /**
     * @brief  Solve the linear system
     *
     * @param  B is the RHS vector
     * @param  X is the solution vector
     * @return void
     */
<<<<<<< HEAD
    template <typename T, typename U>
    void solve(const T& b, U& x) {
=======
    template <typename T, 
      typename U
    >
    void solve(const T& b,
      U& x
    ) {
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
      this->underlying().solve(b, x);
    }


    int getMaxIterations() {
      return this->underlying().getMaxIterations();
    }


    void setMaxIterations(int maxIters) {
      this->underlying().setMaxIterations(maxIters);
    }


    double getTolerance() {
      return this->underlying().getTolerance();
    }
 

    void setTolerance(double tolerance) {
      this->underlying().setTolerance(tolerance);
    }


  protected:

    LinearIterativeSolverBase() = default;
    ~LinearIterativeSolverBase() = default;

<<<<<<< HEAD
=======
    LinearIterativeSolverBase(LinearIterativeSolverBase&&) = delete;
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
    LinearIterativeSolverBase(const LinearIterativeSolverBase&) = delete;


  private:

    Derived& underlying() {
      return static_cast<Derived&>(*this);
    }
  

    Derived const& underlying() const {
      return static_cast<Derived const&>(*this);
    }  
};

} //end namespace solvers

#endif
