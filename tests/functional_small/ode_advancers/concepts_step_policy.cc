
#include <gtest/gtest.h>
#include "pressio/ode_advancers.hpp"

namespace{
  using namespace pressio::ode;

  struct BasicPolicy1{
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSize<double> & /*unused*/) const{}
  };

  struct BasicPolicy2{
    // wrong order
    void operator()(StepCount /*unused*/,
		    StepSize<double> & /*unused*/,
		    StepStartAt<double> /*unused*/) const{}
  };

  struct BasicPolicy3{
    // dt is taken by value
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSize<double> /*unused*/) const{}
  };

  struct BasicPolicy4{
    // dt is rval ref
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSize<double> && /*unused*/) const{}
  };

  struct Policy1{
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSize<double> & /*unused*/,
		    StepSizeMin<double> & /*unused*/,
		    StepSizeReduction<double> & /*unused*/) const{}
  };

  struct Policy2{
    // wrong order
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSizeMin<double> & /*unused*/,
		    StepSize<double> & /*unused*/,
		    StepSizeReduction<double> & /*unused*/) const{}
  };

  struct Policy3{
    // wrong order
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSizeMin<double> & /*unused*/,
		    StepSizeReduction<double> & /*unused*/,
		    StepSize<double> & /*unused*/) const {}
  };

  struct Policy4{
    // missing ref
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSize<double> /*unused*/,
		    StepSizeMin<double> & /*unused*/,
		    StepSizeReduction<double> & /*unused*/) const{}
  };

  struct Policy5{
    // missing ref
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSize<double> & /*unused*/,
		    StepSizeMin<double> /*unused*/,
		    StepSizeReduction<double> & /*unused*/) const{}
  };

  struct Policy6{
    // missing ref
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSize<double> & /*unused*/,
		    StepSizeMin<double> & /*unused*/,
		    StepSizeReduction<double> /*unused*/) const{}
  };

  struct Policy7{
    // missing ref
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSize<double> /*unused*/,
		    StepSizeMin<double> /*unused*/,
		    StepSizeReduction<double> & /*unused*/) const{}
  };

  struct Policy8{
    // missing ref
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    StepSize<double> /*unused*/,
		    StepSizeMin<double> & /*unused*/,
		    StepSizeReduction<double> /*unused*/) const{}
  };

} //end anonym namespace

TEST(ode, concepts_time_step_size_policy)
{
  using namespace pressio::ode;

  static_assert( StepSizePolicy<BasicPolicy1, double>::value, "");
  static_assert(!StepSizePolicy<BasicPolicy2, double>::value, "");
  static_assert(!StepSizePolicy<BasicPolicy3, double>::value, "");
  static_assert(!StepSizePolicy<BasicPolicy4, double>::value, "");
}

TEST(ode, concepts_time_step_size_policy_with_reduction)
{
  using namespace pressio::ode;

  using time_type = double;
  static_assert(  StepSizePolicyWithReductionScheme<Policy1, time_type>::value, "");
  static_assert( !StepSizePolicyWithReductionScheme<Policy2, time_type>::value, "");
  static_assert( !StepSizePolicyWithReductionScheme<Policy3, time_type>::value, "");
  static_assert( !StepSizePolicyWithReductionScheme<Policy4, time_type>::value, "");
  static_assert( !StepSizePolicyWithReductionScheme<Policy5, time_type>::value, "");
  static_assert( !StepSizePolicyWithReductionScheme<Policy6, time_type>::value, "");

  // these cannot be on yet because of the template ambigous error
  //static_assert( !time_step_size_policy_with_step_local_reduction<Policy7, time_type>::value, "");
  //static_assert( !time_step_size_policy_with_step_local_reduction<Policy8, time_type>::value, "");
}
