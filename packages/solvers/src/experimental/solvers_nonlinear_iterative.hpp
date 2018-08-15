
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

  	

  protected:

  	NonLinearIterativeSolver() : base_type() {}

};
	
} // end namespace solvers

#endif