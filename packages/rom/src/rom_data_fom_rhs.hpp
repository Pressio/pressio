
#ifndef ROM_DATA_FOM_RHS_HPP_
#define ROM_DATA_FOM_RHS_HPP_

#include "rom_ConfigDefs.hpp"
#include "rom_fwd.hpp"

namespace pressio{ namespace rom{

template <typename fom_rhs_type>
struct FomRhsData<fom_rhs_type>
{
  using fom_rhs_t = fom_rhs_type;
  FomRhsData() = delete;

  template <
    typename _fom_rhs_type = fom_rhs_type,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_vector_wrapper<_fom_rhs_type>::value
      > * = nullptr
    >
  FomRhsData(const _fom_rhs_type & fomRhsIn)
    : fomRhs_{fomRhsIn}
  {
    ::pressio::containers::ops::set_zero(fomRhs_);
  }

#ifdef HAVE_PYBIND11
  template <
    typename _fom_rhs_type = fom_rhs_type,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_cstyle_array_pybind11<_fom_rhs_type>::value
      > * = nullptr
    >
  FomRhsData(const _fom_rhs_type & fomRhsIn)
    : fomRhs_{{_fom_rhs_type(const_cast<_fom_rhs_type &>(fomRhsIn).request())}}
  {
    ::pressio::containers::ops::set_zero(fomRhs_);
  }
#endif

  ~FomRhsData() = default;

protected:
  mutable fom_rhs_t fomRhs_ = {};
};//end class

}}//end namespace pressio::rom
#endif
