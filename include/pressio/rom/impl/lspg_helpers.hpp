
#ifndef PRESSIO_ROM_IMPL_LSPG_HELPERS_HPP_
#define PRESSIO_ROM_IMPL_LSPG_HELPERS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template<class T= void>
void valid_scheme_for_lspg_else_throw(::pressio::ode::StepScheme name){
  if (   name != ::pressio::ode::StepScheme::BDF1
      && name != ::pressio::ode::StepScheme::BDF2)
  {
    throw std::runtime_error("LSPG currently accepting BDF1 or BDF2");
  }
}

}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_LSPG_HELPERS_HPP_
