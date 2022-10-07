
#include <gtest/gtest.h>
#include "pressio/ode_concepts.hpp"

using state_type = std::vector<float>;

namespace{
  using namespace pressio::ode;

  struct S1{
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    state_type & /*unused*/) const{}
  };

  struct S2{
    // wrong order
    void operator()(StepCount /*unused*/,
		    state_type & /*unused*/,
		    StepStartAt<double> /*unused*/) const{}
  };

  struct S3{
    // wrong order
    void operator()(state_type & /*unused*/,
		    StepCount /*unused*/,
		    StepStartAt<double> /*unused*/) const{}
  };

  struct S4{
    // missing ref
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    state_type /*unused*/) const{}
  };

  struct S5{
    // rval ref
    void operator()(StepCount /*unused*/,
		    StepStartAt<double> /*unused*/,
		    state_type && /*unused*/) const{}
  };

} //end anonym namespace

TEST(ode, concepts_state_observer)
{
  using namespace pressio::ode;

#ifdef PRESSIO_ENABLE_CXX20
  static_assert( StateGuesser<S1, double, state_type>, "");
  static_assert( !StateGuesser<S2, double, state_type>, "");
  static_assert( !StateGuesser<S3, double, state_type>, "");
  static_assert( !StateGuesser<S4, double, state_type>, "");
  static_assert( !StateGuesser<S5, double, state_type>, "");
#else
  static_assert( StateGuesser<S1, double, state_type>::value, "");
  static_assert( !StateGuesser<S2, double, state_type>::value, "");
  static_assert( !StateGuesser<S3, double, state_type>::value, "");
  static_assert( !StateGuesser<S4, double, state_type>::value, "");
  static_assert( !StateGuesser<S5, double, state_type>::value, "");
#endif  
}
