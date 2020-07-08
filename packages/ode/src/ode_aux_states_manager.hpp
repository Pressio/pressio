
#ifndef ODE_AUX_STATES_CONTAINER_HPP_
#define ODE_AUX_STATES_CONTAINER_HPP_

namespace pressio{ namespace ode{

template<typename T, std::size_t n>
class AuxStatesManager;

// partially specialize for implicit scheme
template<typename T, std::size_t n>
class AuxStatesManager
{
public:
  using data_type = ::pressio::containers::StaticCollection<T, n>;

  template <typename ... Args>
  AuxStatesManager(Args && ... args)
    : data_( std::forward<Args>(args)... ){}

  ~AuxStatesManager() = default;

public:
  static constexpr std::size_t size() {
    return data_type::size();
  }

  // n-1
  T & get(ode::nMinusOne tag){ return data_(0); }
  T const & get(ode::nMinusOne tag) const{ return data_(0); }

  // n-2
  T & get(ode::nMinusTwo tag){ return data_(1); }
  T const & get(ode::nMinusTwo tag) const{ return data_(1); }

  // n-3
  T & get(ode::nMinusThree tag){ return data_(2); }
  T const & get(ode::nMinusThree tag) const{ return data_(2); }

  // n-4
  T & get(ode::nMinusFour tag){ return data_(3); }
  T const & get(ode::nMinusFour tag) const{ return data_(3); }

private:
  data_type data_;
};

}}//end namespace pressio::ode
#endif
