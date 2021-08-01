
#include <gtest/gtest.h>
#include "pressio_ode_implicit.hpp"

TEST(ode_implicit, staticCheckDtSetter){
  using namespace pressio;
  using step_t = ode::step_type;
  using time_t = double;

  {
    const auto lambda = [](const step_t &, const time_t &, time_t &){};
    using l_t = decltype(lambda);
    static_assert( ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, time_t &) -> void{};
    using l_t = decltype(lambda);
    static_assert( ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, time_t &) -> double{ return 11.0; };
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, time_t &, time_t & ){};
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](step_t &, time_t &, time_t & ){};
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, time_t &, const time_t & ){};
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](time_t &, const time_t & ){};
    using l_t = decltype(lambda);
    static_assert( !ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

}



TEST(ode_implicit, staticCheckGuesserFunctor){
  using namespace pressio;
  using step_t = ode::step_type;
  using time_t = double;

  {
    struct Guesser{
      void operator()(const step_t &, const time_t &, time_t &) const
      {}
    };
    using l_t = Guesser;;
    static_assert( ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    struct Guesser{
      double operator()(const step_t &, const time_t &, time_t &)
      {return 0.0;}
    };
    using l_t = Guesser;;
    static_assert( !ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(const step_t &, time_t &, time_t & )
      {}
    };
    using l_t = Guesser;;
    static_assert( !ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(step_t &, time_t &, time_t & )
      {}
    };
    using l_t = Guesser;;
    static_assert( !ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(const step_t &, time_t &, const time_t & )
      {}
    };
    using l_t = Guesser;;
    static_assert( !ode::constraints::time_step_size_manager<l_t, step_t, time_t>::value, "" );
  }

}
