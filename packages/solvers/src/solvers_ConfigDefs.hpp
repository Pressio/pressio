
#ifndef SOLVERS_CONFIGDEFS_HPP_
#define SOLVERS_CONFIGDEFS_HPP_

#include "solvers_config.h"
#include "../../core/src/core_ConfigDefs.hpp"
#include <iostream>
#include <cassert>

namespace rompp{ namespace solvers{

struct L2Norm{};
struct L1Norm{};

namespace iterative{


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


namespace gn{

struct noLineSearch{};
struct ArmijoLineSearch{};

}//end namespace rompp::solvers::gn



}}}//end namespace rompp::solvers::iterative
#endif
