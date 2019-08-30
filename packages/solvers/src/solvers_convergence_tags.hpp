
#ifndef SOLVERS_CONVERGENCE_TAGS_HPP_
#define SOLVERS_CONVERGENCE_TAGS_HPP_

#include "solvers_norm_tags.hpp"

namespace pressio{ namespace solvers{ namespace iterative{

namespace converged_when{

template <typename norm_type>
struct absoluteNormCorrectionBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for absoluteNormCorrectionBelowTol, it must be L1Norm or L2Norm");
};

template <typename norm_type>
struct absoluteNormResidualBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for absoluteNormResidualBelowTol, it must be L1Norm or L2Norm");
};

template <typename norm_type>
struct relativeNormResidualBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for relativeNormResidualBelowTol, it must be L1Norm or L2Norm");
};

template <typename norm_type>
struct absoluteNormProjectedResidualBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for absoluteNormProjectedResidualBelowTol, it must be L1Norm or L2Norm");
};

template <typename norm_type>
struct relativeNormProjectedResidualBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for relativeNormProjectedResidualBelowTol, it must be L1Norm or L2Norm");
};


struct completingNumMaxIters{};

}//end namespace pressio::solvers::convergedWhen

using default_convergence
	= converged_when::absoluteNormCorrectionBelowTol<L2Norm>;

}}}//end namespace pressio::solvers::iterative
#endif
