#ifndef SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP

struct SolversNonLinearIterativeLeastSquareLevenbergMarquardtPolicy {

  template <
    typename SolverT,
    typename PrecT,
    typename NormT,
    typename SystemT,
    typename VectorT,
    typename std::enable_if<
      core::details::traits<VectorT>::is_vector &&
      solvers::meta::are_vector_compatible<
        typename details::system_traits<SystemT>::vector_type,
        VectorT
      >::value,
      int
    >::type* = nullptr
  >
  static auto solve(
    const SystemT& sys,
    const VectorT& x0,
    core::defaultTypes::uint maxIterations,
    core::defaultTypes::uint maxNonLinearIterations,
    double tolerance,
    double nonLinearTolerance
  ) {

    auto dy = sys.residual(x0);
    auto Ja = sys.jacobian(x0);

    auto solver = LinearSolvers::createIterativeSolver<SolverT, typename SystemT::matrix_type, PrecT>();
    solver.setMaxIterations(maxIterations);
    solver.setTolerance(tolerance);

    auto xOld = x0;
    auto xNew = x0;
    core::defaultTypes::uint iStep = 1;

    do {
      auto A =  

      solver.resetLinearSystem(Ja);
      xNew = xNew - solver.solve(dy);
    } while (iStep++ < maxNonLinearIterations && NormT::template compute_norm_difference(xOld, xNew) > nonLinearTolerance);


    return xNew;
  }

};

#endif
