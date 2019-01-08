
#ifndef ROM_DATA_BASE_HPP_
#define ROM_DATA_BASE_HPP_

#include "rom_ConfigDefs.hpp"

namespace rompp{ namespace rom{


template <typename app_res_w_type,
	  int maxNrhs>
class RomRHSData{

protected:
  // default constructor
  RomRHSData() = delete;

  // non-default constructors
  template <int _maxNrhs = maxNrhs,
            core::meta::enable_if_t<_maxNrhs==1> * = nullptr>
  RomRHSData(const app_res_w_type & r0fom)
    : appRHS_{r0fom}{}

  ~RomRHSData() = default;

protected:
  mutable std::array<app_res_w_type, maxNrhs> appRHS_ = {};

};//end class





template <typename app_state_w_type,
	  typename phi_op_type,
	  int maxNstates>
class RomStateData{

protected:

  // default constructor
  RomStateData() = delete;


  // non-default constructors
  template <int _maxNstates = maxNstates,
            core::meta::enable_if_t<_maxNstates==0> * = nullptr>
  RomStateData(const app_state_w_type & y0fom, phi_op_type & phiOp)
    : y0FOM_(&y0fom),
      yFOM_(y0fom),
      phi_(&phiOp){}

  // non-default constructors
  template <int _maxNstates = maxNstates,
            core::meta::enable_if_t<_maxNstates==1> * = nullptr>
  RomStateData(const app_state_w_type & y0fom, phi_op_type & phiOp)
    : y0FOM_(&y0fom),
      yFOM_(y0fom),
      yFOMold_{y0fom},
      phi_(&phiOp){}

  // non-default constructors
  template <int _maxNstates = maxNstates,
            core::meta::enable_if_t<_maxNstates==2> * = nullptr>
  RomStateData(const app_state_w_type & y0fom, phi_op_type & phiOp)
    : y0FOM_(&y0fom),
      yFOM_(y0fom),
      yFOMold_{y0fom, y0fom},
      phi_(&phiOp){}

  ~RomStateData() = default;

protected:

  template <typename ode_state_t>
  void reconstructCurrentFOMState(const ode_state_t & odeY) const{
    phi_->apply(odeY, yFOM_);
    yFOM_ += (*y0FOM_);
  }
  //----------------------------------------------------------------


  template <int numAuxStates,
	    typename ode_state_t,
            core::meta::enable_if_t<numAuxStates==1> * = nullptr>
  void reconstructFOMStates(
	  const ode_state_t & odeY,
          const std::array<ode_state_t,numAuxStates> & odeYprev) const{

    reconstructCurrentFOMState(odeY);
    phi_->apply(odeYprev[0], yFOMold_[0]);
    yFOMold_[0] += (*y0FOM_);
  }
  //----------------------------------------------------------------

  template <int numAuxStates,
	    typename ode_state_t,
            core::meta::enable_if_t<numAuxStates==2> * = nullptr>
  void reconstructFOMStates(
	  const ode_state_t & odeY,
          const std::array<ode_state_t,numAuxStates> & odeYprev) const{

    reconstructCurrentFOMState(odeY);

    phi_->apply(odeYprev[0], yFOMold_[0]);
    phi_->apply(odeYprev[1], yFOMold_[1]);
    yFOMold_[0] += (*y0FOM_);
    yFOMold_[1] += (*y0FOM_);
  }
  //----------------------------------------------------------------

protected:
  const app_state_w_type * y0FOM_ 			     = nullptr;
  mutable app_state_w_type yFOM_                             = {};
  mutable std::array<app_state_w_type, maxNstates> yFOMold_  = {};
  phi_op_type * phi_                                         = nullptr;

};//end class


}}//end namespace rompp::rom
#endif
