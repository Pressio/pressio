
#ifndef SOLVERS_NONLINEAR_BASE_HPP
#define SOLVERS_NONLINEAR_BASE_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../meta/solvers_is_legitimate_system_for_nonlinear_solver.hpp"
#include "../../../containers/src/vector/containers_vector_meta.hpp"

namespace pressio{ namespace solvers{

/**
 * @brief Base class for nonlinear solver
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a nonlinear solver.
 */
template <typename Derived>
struct NonLinearSolverBase {

  NonLinearSolverBase()	 = default;
  ~NonLinearSolverBase() = default;
  NonLinearSolverBase(const NonLinearSolverBase &) = delete;

  std::string getConvergenceConditionDescription() const{
    return convergenceConditionDescription_;
  }

  template <typename system_t, typename state_t >
  void solve(const system_t & sys, state_t & x){
    // static_assert( ::pressio::solvers::meta::is_legitimate_system_for_nonlinear_solver<
    //   system_t>::value,
    // 		   "The system obj type is not valid for non-linear solver");

    static_cast<Derived&>(*this).solveImpl(sys, x);
  }

private:
  std::string convergenceConditionDescription_ = "none";

  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<Derived>::type;
};

}}//end namespace pressio::solvers
#endif
