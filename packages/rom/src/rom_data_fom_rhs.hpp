
#ifndef ROM_FOM_RHS_DATA_HPP_
#define ROM_FOM_RHS_DATA_HPP_

#include "rom_ConfigDefs.hpp"

namespace rompp{ namespace rom{

template <typename fom_rhs_w_type>
struct FomRhsData{

  using fom_rhs_w_t = fom_rhs_w_type;

  FomRhsData() = delete;

  FomRhsData(const fom_rhs_w_t & fomRhs0)
    : fomRhs_{fomRhs0}{}

  ~FomRhsData() = default;

protected:
  mutable fom_rhs_w_t fomRhs_ = {};

};//end class

}}//end namespace rompp::rom
#endif
