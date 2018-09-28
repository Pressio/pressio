
#ifndef SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP

#include <iostream>
#include <type_traits>

#include "solvers_system_traits.hpp"
#include "solvers_nonlinear_traits.hpp"
#include "solvers_meta_static_checks.hpp"
#include "../solvers_ConfigDefs.hpp"
#include "../../../core/src/meta/core_meta_detection_idiom.hpp"


namespace rompp{
namespace solvers{

namespace todeprecate{

/**
 * Return the transpose of a matrix. To deprecate and implement in core
 **/
template <
  typename MatrixT,
  typename core::meta::enable_if_t<
    core::details::traits<MatrixT>::wrapped_package_identifier == core::details::WrappedPackageIdentifier::Eigen
  >* = nullptr
>
auto transpose(const MatrixT& A)
{
  MatrixT res(A.data()->transpose());
  return res;
}

} // end namespace todeprecate


struct SolversNonLinearIterativeLeastSquareLevenbergMarquardtPolicy {

  template <
    typename SolverT,
    typename PrecT,
    typename NormT,
    typename SystemT,
    typename VectorT,
    typename core::meta::enable_if_t<
      core::details::traits<VectorT>::is_vector &&
      solvers::meta::are_vector_compatible<
        typename details::system_traits<SystemT>::vector_type,
        VectorT
      >::value
    >* = nullptr
  >
  static auto solve(
    const SystemT& sys,
    const VectorT& x0,
    core::defaultTypes::uint maxIterations,
    core::defaultTypes::uint maxNonLinearIterations,
    double tolerance,
    double nonLinearTolerance,
    double lambda
  ) {

    auto dy = sys.residual(x0);
    auto Ja = sys.jacobian(x0);
    auto JaT = todeprecate::transpose(Ja);
    auto lId = decltype(Ja)(Ja.cols(), Ja.cols()).setIdentity();
    lId.addToDiagonal(lambda);

    auto b = JaT * dy;
    auto A = JaT * Ja + lId;

    auto solver = LinearSolvers::createIterativeSolver<SolverT, typename SystemT::matrix_type, PrecT>(A);
    solver.setMaxIterations(maxIterations);
    solver.setTolerance(tolerance);

    auto xOld = x0;
    auto xNew = x0;

    double normO = 0.0;
    double normN = 0.0;
    core::defaultTypes::uint iStep = 1;

    while (iStep++ < maxNonLinearIterations) {
      xNew = xOld - solver.solve(b);
      auto dyN = sys.residual(xNew);

      normO = NormT::template compute_norm(dy);
      normN = NormT::template compute_norm(dyN);

      if (normN >= normO) {
        // Step not accepted
        lambda *= 2;
        lId.setIdentity();
        lId.addToDiagonal(lambda);

        A = JaT * Ja + lId;
        solver.resetLinearSystem(A);
      } else {
        // Step accepted
        if ((normO-normN) < nonLinearTolerance) {break;}
        lambda *= 0.8;

        xOld = xNew;
        dy = sys.residual(xOld);
        Ja = sys.jacobian(xOld);
        JaT = todeprecate::transpose(Ja);
        lId.setIdentity();
        lId.addToDiagonal(lambda);

        b = JaT * dy;
        A = JaT * Ja + lId;
        solver.resetLinearSystem(A);
      }
    }
    return xNew;
  }

};

} // end namespace solvers
} // end namespace rompp

#endif
