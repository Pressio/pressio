
#ifndef ROM_DATA_FOM_RHS_HPP_
#define ROM_DATA_FOM_RHS_HPP_

#include "rom_ConfigDefs.hpp"

namespace rompp{ namespace rom{


template <
  typename fom_rhs_type,
  typename enable = void>
struct FomRhsData;


template <typename fom_rhs_type>
struct FomRhsData<
  fom_rhs_type,
  mpl::enable_if_t<
    ::rompp::algebra::meta::is_algebra_vector_wrapper<fom_rhs_type>::value
    >
  >
{
  using fom_rhs_t = fom_rhs_type;
  FomRhsData() = delete;

  FomRhsData(const fom_rhs_type & fomRhs0)
    : fomRhs_{fomRhs0}{
    // reset to zero to be safe
    fomRhs_.setZero();
  }

  ~FomRhsData() = default;

protected:
  mutable fom_rhs_t fomRhs_ = {};
};//end class





#ifdef HAVE_PYBIND11
template <typename fom_rhs_type>
struct FomRhsData<
  fom_rhs_type,
  mpl::enable_if_t<
    ::rompp::algebra::meta::is_cstyle_array_pybind11<fom_rhs_type>::value
    >
  >
{
  using fom_rhs_t = fom_rhs_type;
  FomRhsData() = delete;

  FomRhsData(const fom_rhs_type & fomRhs0)
    : fomRhs_{{fom_rhs_type(const_cast<fom_rhs_type &>(fomRhs0).request())}}
  {
    // reset to zero to be safe
    ::rompp::algebra::ops::set_zero(fomRhs_);
    std::cout << "FomRhsData:: fomRhs_" << fomRhs_.data() << std::endl;
  }

  ~FomRhsData() = default;

protected:
  mutable fom_rhs_t fomRhs_ = {};
};//end class
#endif


}}//end namespace rompp::rom
#endif
