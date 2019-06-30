
#ifndef SOLVERS_IS_NON_DEFAULT_CONVERGENCE_TAG_HPP_
#define SOLVERS_IS_NON_DEFAULT_CONVERGENCE_TAG_HPP_

#include "../solvers_convergence_tags.hpp"
#include "../solvers_norm_tags.hpp"

namespace rompp{ namespace solvers{ namespace meta {

template <typename T, typename enable = void>
struct is_non_default_convergence_tag
  : std::false_type{};

template <typename T>
struct is_non_default_convergence_tag<
  T,
  ::rompp::mpl::enable_if_t<
    std::is_same<
      T,
      iterative::converged_when::absoluteNormCorrectionBelowTol<L2Norm>>::value or
    std::is_same<
      T,
      iterative::converged_when::absoluteNormCorrectionBelowTol<L1Norm>>::value or
    std::is_same<
      T,
      iterative::converged_when::absoluteNormResidualBelowTol<L2Norm>>::value or
    std::is_same<
      T,
      iterative::converged_when::absoluteNormResidualBelowTol<L1Norm>>::value or
    std::is_same<
      T,
      iterative::converged_when::relativeNormResidualBelowTol<L2Norm>>::value or
    std::is_same<
      T,
    iterative::converged_when::relativeNormResidualBelowTol<L1Norm>>::value or
    std::is_same<
    T,
    iterative::converged_when::completingNumMaxIters>::value
    >
  > : std::true_type{};

}}} // namespace rompp::solvers::meta
#endif
