#ifndef ROM_UTILS_SET_ROM_COEFFICIENTS_HPP_
#define ROM_UTILS_SET_ROM_COEFFICIENTS_HPP_
namespace pressio{ namespace rom{ namespace utils{ 

template <typename linear_solver_t,
          typename basis_t, 
          typename rom_state_t, 
          typename fom_state_t>
void setRomCoefficientsL2Projection(const basis_t & phi,
                      rom_state_t & romState,
                      const fom_state_t & yFOM_IC,
                      const fom_state_t & yRef, 
                      const int & romSize)
{
  /* 
  Function to compute the ROM coefficients from optimal 
  L^2 projection of yFOM
  */


  using hessian_t = typename linear_solver_t::matrix_type;
  //initialize linear solver. This is only used here
  linear_solver_t linear_solver;
  // create the system matrix, phi^T phi
  hessian_t H(romSize,romSize);
  ::pressio::containers::ops::dot(phi,phi,H);
  //create a vector to store yFOM - yRef
  fom_state_t b(yFOM_IC);
  pressio::containers::ops::do_update(b,::pressio::utils::constants::one(),yRef,::pressio::utils::constants::negOne());
  // compute phi^T b
  rom_state_t r(romSize);
  pressio::containers::ops::dot(phi, b, r);
  //solve system for optimal L2 projection
  linear_solver.solveAllowMatOverwrite(H, r, romState);
}

}}} //end namespace pressio

#endif
