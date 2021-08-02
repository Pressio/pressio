
#include <gtest/gtest.h>
#include "pressio_ode.hpp"

namespace
{

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

//====================================================================

struct ValidSystemWithVelocity1{
  using scalar_type = float;
  using state_type = std::vector<float>;
  using velocity_type = state_type;
  velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
};

struct InvalidSystemWithVelocity1{
  using scalar_type = float;
  using state_type = std::vector<float>;
  using velocity_type = state_type;
  // velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
};

struct InvalidSystemWithVelocity2{
  // using scalar_type = float;
  using state_type = std::vector<float>;
  using velocity_type = state_type;
  velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
};

struct ValidSystemWithVelocityAndJacobian1{
  using scalar_type = float;
  using state_type = std::vector<float>;
  using velocity_type = state_type;
  using jacobian_type = std::vector<std::vector<float>>;
  velocity_type createVelocity() const{ return velocity_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
  void jacobian(const state_type &, double time, jacobian_type &) const{}
};

struct ValidDiscreteTimeSystem{
  using scalar_type = float;
  using state_type = std::vector<float>;
  using discrete_time_residual_type = state_type;
  using discrete_time_jacobian_type = std::vector<std::vector<float>>;

  discrete_time_residual_type createDiscreteTimeResidual() const{ 
    return discrete_time_residual_type(); }
  discrete_time_jacobian_type createDiscreteTimeJacobian() const{ 
    return discrete_time_jacobian_type(); }

  void discreteTimeResidual(pressio::ode::step_count_type, 
                              double time, 
                              double dt, 
                              discrete_time_residual_type &, 
                              const state_type &) const{}
  void discreteTimeResidual(pressio::ode::step_count_type, 
                              double time, 
                              double dt, 
                              discrete_time_residual_type &, 
                              const state_type &,
                              const state_type &) const{}
  void discreteTimeJacobian(pressio::ode::step_count_type, 
                              double time, 
                              double dt, 
                              discrete_time_jacobian_type &, 
                              const state_type &) const{}
  void discreteTimeJacobian(pressio::ode::step_count_type, 
                              double time, 
                              double dt, 
                              discrete_time_jacobian_type &, 
                              const state_type &,
                              const state_type &) const{}
};

} //end namespace

TEST(ode, concepts_collector)
{
  using state_type = std::vector<double>;
  using time_type = double;
  static_assert(pressio::ode::collector<ValidCollector1<state_type>, time_type, state_type>::value, "");
  static_assert(pressio::ode::collector<ValidCollector2<state_type>, time_type, state_type>::value, "");
  static_assert(pressio::ode::collector<ValidCollector3<state_type>, time_type, state_type>::value, "");
  static_assert(pressio::ode::collector<ValidCollector4<state_type>, time_type, state_type>::value, "");
  static_assert(!pressio::ode::collector<InvalidCollector1<state_type>, time_type, state_type>::value, "");
  static_assert(!pressio::ode::collector<InvalidCollector2<state_type>, time_type, state_type>::value, "");

  auto myc = [](pressio::ode::step_count_type, double time, const state_type &){};
  static_assert(pressio::ode::collector<decltype(myc), time_type, state_type>::value, "");
}

TEST(ode, concepts_continuous_time_system)
{
  using namespace pressio::ode;
  static_assert(continuous_time_system_with_at_least_velocity<ValidSystemWithVelocity1>::value, "");
  static_assert(!continuous_time_system_with_at_least_velocity<InvalidSystemWithVelocity1>::value, "");
  static_assert(!continuous_time_system_with_at_least_velocity<InvalidSystemWithVelocity2>::value, "");
  static_assert(continuous_time_system_with_at_least_velocity<ValidSystemWithVelocityAndJacobian1>::value, "");

  static_assert(continuous_time_system_with_user_provided_jacobian<ValidSystemWithVelocityAndJacobian1>::value, "");
  static_assert(!continuous_time_system_with_user_provided_jacobian<ValidSystemWithVelocity1>::value, "");
  static_assert(!continuous_time_system_with_user_provided_jacobian<ValidDiscreteTimeSystem>::value, "");
}

TEST(ode, concepts_discrete_time_system)
{
  using namespace pressio::ode;
  static_assert(discrete_time_system_with_user_provided_jacobian<ValidDiscreteTimeSystem>::value, "");
  static_assert(!discrete_time_system_with_user_provided_jacobian<ValidSystemWithVelocity1>::value, "");
  static_assert(!discrete_time_system_with_user_provided_jacobian<ValidSystemWithVelocityAndJacobian1>::value, "");
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
    static_assert( ode::time_step_size_manager<decltype(lambda), step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, time_t &) -> void{};
    static_assert( ode::time_step_size_manager<decltype(lambda), step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, const time_t &, time_t &) -> double{ return 11.0; };
    static_assert( !ode::time_step_size_manager<decltype(lambda), step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, time_t &, time_t & ){};
    static_assert( !ode::time_step_size_manager<decltype(lambda), step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](step_t &, time_t &, time_t & ){};
    static_assert( !ode::time_step_size_manager<decltype(lambda), step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](const step_t &, time_t &, const time_t & ){};
    static_assert( !ode::time_step_size_manager<decltype(lambda), step_t, time_t>::value, "" );
  }

  {
    const auto lambda = [](time_t &, const time_t & ){};
    static_assert( !ode::time_step_size_manager<decltype(lambda), step_t, time_t>::value, "" );
  }

  {
    struct Setter{
      void operator()(const step_t &, const time_t &, time_t &) const{}
    };    
    static_assert( ode::time_step_size_manager<Setter, step_t, time_t>::value, "" );
  }

  {
    struct Setter{
      double operator()(const step_t &, const time_t &, time_t &){return 0.0;}
    };    
    static_assert( !ode::time_step_size_manager<Setter, step_t, time_t>::value, "" );
  }

  {
    struct Setter{
      void operator()(const step_t &, time_t &, time_t & ){}
    };    
    static_assert( !ode::time_step_size_manager<Setter, step_t, time_t>::value, "" );
  }

  {
    struct Setter{
      void operator()(step_t &, time_t &, time_t & ){}
    };    
    static_assert( !ode::time_step_size_manager<Setter, step_t, time_t>::value, "" );
  }

  {
    struct Setter{
      void operator()(const step_t &, time_t &, const time_t & ){}
    };    
    static_assert( !ode::time_step_size_manager<Setter, step_t, time_t>::value, "" );
  }
}
