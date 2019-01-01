
#ifndef SOLVERS_ITERATIVE_BASE_HPP
#define SOLVERS_ITERATIVE_BASE_HPP

#include "../solvers_ConfigDefs.hpp"

namespace rompp { namespace solvers {

template<typename SolverT, typename scalar_t>
struct IterativeBase
{
  using uint_t = core::default_types::uint;

  inline uint_t getMaxIterations() {
    return maxIters_;
  }

  void setMaxIterations(uint_t maxIters) {
    maxIters_ = maxIters;
  }

  inline scalar_t getTolerance() {
    return tolerance_;
  }

  void setTolerance(scalar_t tolerance) {
    tolerance_ = tolerance;
  }

public:

  IterativeBase() = default;
  IterativeBase(const IterativeBase &) = delete;
  ~IterativeBase() = default;

  // IterativeBase(std::shared_ptr<SolverT> solver) :
  //   base_type(solver), maxIters_(100), tolerance_(1.0e-6) {};

protected:
  uint_t maxIters_    = static_cast<uint_t>(100);
  scalar_t tolerance_ = static_cast<scalar_t>(0.000001);
};


}} //end namespace rompp::solvers

#endif
