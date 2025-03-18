
#include <gtest/gtest.h>
#include "pressio/ode_concepts.hpp"

using state_type = std::vector<float>;

namespace{
  using namespace pressio::ode;

  struct Obs1{
    void operator()(StepCount /*unused*/,
		    double /*unused*/,
		    const state_type & /*unused*/) const{}
  };

  struct Obs2{
    // wrong order
    void operator()(StepCount /*unused*/,
		    const state_type & /*unused*/,
		    double /*unused*/) const{}
  };

  struct Obs3{
    // wrong order
    void operator()(const state_type & /*unused*/,
		    StepCount /*unused*/,
		    double /*unused*/) const{}
  };

} //end anonym namespace

TEST(ode, concepts_state_observer)
{
  using namespace pressio::ode;

  static_assert( StateObserver<Obs1, double, state_type>::value, "");
  static_assert( !StateObserver<Obs2, double, state_type>::value, "");
  static_assert( !StateObserver<Obs3, double, state_type>::value, "");

}
