
#ifndef SOLVERS_CONVERGENCE_TAGS_HPP_
#define SOLVERS_CONVERGENCE_TAGS_HPP_

#include "solvers_norm_tags.hpp"

namespace rompp{ namespace solvers{ namespace iterative{

namespace converged_when{

template <typename norm_type>
struct absoluteNormCorrectionBelowTol{
  using norm_t = norm_type;
};

template <typename norm_type>
struct absoluteNormResidualBelowTol{
  using norm_t = norm_type;
};

template <typename norm_type>
struct relativeNormResidualBelowTol{
  using norm_t = norm_type;
};

struct completingNumMaxIters{};

}//end namespace rompp::solvers::convergedWhen

using default_convergence
	= converged_when::absoluteNormCorrectionBelowTol<L2Norm>;

}}}//end namespace rompp::solvers::iterative
#endif
