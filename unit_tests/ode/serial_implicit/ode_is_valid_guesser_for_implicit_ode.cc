
#include <gtest/gtest.h>
#include "pressio_ode_implicit.hpp"

TEST(ode_implicit, staticCheckGuesserLambda){
  using namespace pressio;
  using sc_t = double;
  using step_t = ode::types::step_t;
  using time_t = sc_t;
  using state_t = std::vector<sc_t>;

  {
    const auto lambda = [](const step_t &, const time_t &, state_t &){};
    using l_t = decltype(lambda);
    static_assert( ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, state_t &) -> void{};
    using l_t = decltype(lambda);
    static_assert( ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, state_t &) -> double{ return 11.0; };
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, time_t &, state_t & ){};
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](step_t &, time_t &, state_t & ){};
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, state_t &, const time_t & ){};
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }


  {
    const auto lambda = [](state_t &, const time_t & ){};
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

}



TEST(ode_implicit, staticCheckGuesserFunctor){
  using namespace pressio;
  using sc_t = double;
  using step_t = ode::types::step_t;
  using time_t = sc_t;
  using state_t = std::vector<sc_t>;

  {
    struct Guesser{
      void operator()(const step_t &, const time_t &, state_t &) const
      {}
    };
    using l_t = Guesser;;
    static_assert( ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

  {
    struct Guesser{
      double operator()(const step_t &, const time_t &, state_t &)
      {return 0.0;}
    };
    using l_t = Guesser;;
    static_assert( !ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(const step_t &, time_t &, state_t & )
      {}
    };
    using l_t = Guesser;;
    static_assert( !ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(step_t &, time_t &, state_t & )
      {}
    };
    using l_t = Guesser;;
    static_assert( !ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(const step_t &, state_t &, const time_t & )
      {}
    };
    using l_t = Guesser;;
    static_assert( !ode::constraints::is_legitimate_guesser<l_t, step_t, time_t, state_t>::value, "" );
  }

}
