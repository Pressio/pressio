
#include <gtest/gtest.h>
#include "pressio/ode_advancers.hpp"

namespace{

template<typename state_type>
struct ValidCollector1{
  void operator()(pressio::ode::step_count_type, double time, const state_type &){}
};

template<typename state_type>
struct ValidCollector2{
  void operator()(pressio::ode::step_count_type, const state_type &, double time){}
};

template<typename state_type>
struct ValidCollector3{
  void operator()(const state_type &, pressio::ode::step_count_type, double time){}
};

template<typename state_type>
struct ValidCollector4{
  void operator()(double time, const state_type &, pressio::ode::step_count_type){}
};

template<typename state_type>
struct InvalidCollector1{
  void operator()(pressio::ode::step_count_type, const state_type &){}
};

template<typename state_type>
struct InvalidCollector2{
  void operator()(const state_type &){}
};

struct SteppableClass
{
  void doStep(const std::vector<float> state,
              const double time, 
              const double dt, 
              const int32_t step)
  {}
};

} //end anonym namespace


TEST(ode, concepts_steppable)
{
  using namespace pressio::ode;
  using state_type = std::vector<float>;
  static_assert(steppable_with<void, SteppableClass, state_type, double>::value, "");
}

TEST(ode, concepts_collector)
{
  using state_type = std::vector<double>;
  using time_type = double;
  static_assert(pressio::ode::observer<ValidCollector1<state_type>, time_type, state_type>::value, "");
  static_assert(pressio::ode::observer<ValidCollector2<state_type>, time_type, state_type>::value, "");
  static_assert(pressio::ode::observer<ValidCollector3<state_type>, time_type, state_type>::value, "");
  static_assert(pressio::ode::observer<ValidCollector4<state_type>, time_type, state_type>::value, "");
  static_assert(!pressio::ode::observer<InvalidCollector1<state_type>, time_type, state_type>::value, "");
  static_assert(!pressio::ode::observer<InvalidCollector2<state_type>, time_type, state_type>::value, "");

  auto myc = [](pressio::ode::step_count_type, double time, const state_type &){};
  static_assert(pressio::ode::observer<decltype(myc), time_type, state_type>::value, "");
}

TEST(ode, concepts_guesser)
{
  using namespace pressio;
  using step_t = pressio::ode::step_count_type;
  using time_t = double;
  using state_t = std::vector<double>;

  {
    const auto lambda = [](const step_t &, const time_t &, state_t ){};    
    static_assert( pressio::ode::is_legitimate_guesser<decltype(lambda), step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, state_t &) -> void{};    
    static_assert( pressio::ode::is_legitimate_guesser<decltype(lambda), step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, state_t &) -> double{ return 11.0; };    
    static_assert( !pressio::ode::is_legitimate_guesser<decltype(lambda), step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, time_t &, state_t & ){};    
    static_assert( !pressio::ode::is_legitimate_guesser<decltype(lambda), step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](step_t &, time_t &, state_t & ){};    
    static_assert( !pressio::ode::is_legitimate_guesser<decltype(lambda), step_t, time_t, state_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, state_t &, const time_t & ){};    
    static_assert( !pressio::ode::is_legitimate_guesser<decltype(lambda), step_t, time_t, state_t>::value, "" );
  }


  {
    const auto lambda = [](state_t &, const time_t &){};    
    static_assert( !pressio::ode::is_legitimate_guesser<decltype(lambda), step_t, time_t, state_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(const step_t &, const time_t &, state_t &) const{}
    };    
    static_assert( pressio::ode::is_legitimate_guesser<Guesser, step_t, time_t, state_t>::value, "" );
  }

  {
    struct Guesser{
      double operator()(const step_t &, const time_t &, state_t &){return 0.0;}
    };
    static_assert( !pressio::ode::is_legitimate_guesser<Guesser, step_t, time_t, state_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(const step_t &, time_t &, state_t & ){}
    };
    static_assert( !pressio::ode::is_legitimate_guesser<Guesser, step_t, time_t, state_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(step_t &, time_t &, state_t & ){}
    };
    static_assert( !pressio::ode::is_legitimate_guesser<Guesser, step_t, time_t, state_t>::value, "" );
  }

  {
    struct Guesser{
      void operator()(const step_t &, state_t &, const time_t & ){}
    };
    static_assert( !pressio::ode::is_legitimate_guesser<Guesser, step_t, time_t, state_t>::value, "" );
  }
}


TEST(ode, concepts_time_step_setter)
{
  using namespace pressio;
  using step_t = pressio::ode::step_count_type;
  using time_t = double;

  {
    const auto lambda = [](const step_t &, const time_t &, time_t &){};
    static_assert( ode::time_step_size_manager<decltype(lambda), time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, time_t &) -> void{};
    static_assert( ode::time_step_size_manager<decltype(lambda), time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, time_t &) -> double{ return 11.0; };
    static_assert( !ode::time_step_size_manager<decltype(lambda), time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, time_t &, time_t & ){};
    static_assert( !ode::time_step_size_manager<decltype(lambda), time_t>::value, "" );
  }

  {
    const auto lambda = [](step_t &, time_t &, time_t & ){};
    static_assert( !ode::time_step_size_manager<decltype(lambda), time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, time_t &, const time_t & ){};
    static_assert( !ode::time_step_size_manager<decltype(lambda), time_t>::value, "" );
  }

  {
    const auto lambda = [](time_t &, const time_t & ){};
    static_assert( !ode::time_step_size_manager<decltype(lambda), time_t>::value, "" );
  }

  {
    struct Setter{
      void operator()(const step_t &, const time_t &, time_t &) const{}
    };    
    static_assert( ode::time_step_size_manager<Setter, time_t>::value, "" );
  }

  {
    struct Setter{
      double operator()(const step_t &, const time_t &, time_t &){return 0.0;}
    };    
    static_assert( !ode::time_step_size_manager<Setter, time_t>::value, "" );
  }

  {
    struct Setter{
      void operator()(const step_t &, time_t &, time_t & ){}
    };    
    static_assert( !ode::time_step_size_manager<Setter, time_t>::value, "" );
  }

  {
    struct Setter{
      void operator()(step_t &, time_t &, time_t & ){}
    };    
    static_assert( !ode::time_step_size_manager<Setter, time_t>::value, "" );
  }

  {
    struct Setter{
      void operator()(const step_t &, time_t &, const time_t & ){}
    };    
    static_assert( !ode::time_step_size_manager<Setter, time_t>::value, "" );
  }
}
