
#include "apps_ks1d_eigen.hpp"

namespace pressio{ namespace apps{ 

void KS1dEigen::velocity(const state_type & u,          
            const scalar_type /* t */,
	velocity_type & rhs) const{

	scalar_type dxInv2 = dxInv_ * dxInv_;

	for (ui_t i=0; i<Nnode_; ++i){

		// Boundary conditions
		scalar_type up = (i < Nnode_ - 1) ? u(i+1) : 0;//u[0];
		scalar_type upp = (i < Nnode_ - 2) ? u(i+2) : u(i);//u[0];
		scalar_type um = (i > 0) ? u(i-1) : 0; //u[N_GRID-1];
		scalar_type umm = (i > 1) ? u(i-2) : u(i);//u[N_GRID-1];

		// du/dx
		scalar_type ux = 0.5 * dxInv_ * (up - um);

		// du^2 / dx
		scalar_type u2x = 0.5 * dxInv_ * (up*up - um*um);

		//d2u / dx2
		scalar_type uxx =  dxInv2 * (up + um - 2 * u(i)) ;

		//d4u / dx4
		scalar_type uxxp = dxInv2 * (upp + u(i) - 2 * up);
        scalar_type uxxm = dxInv2 * (umm + u(i) - 2 * um);
        scalar_type uxxxx = dxInv2 * (uxxp + uxxm - 2 * uxx);

		rhs(i) = -1.0 * mu_(0) * ux - 0.5 * u2x - uxx - mu_(1) * uxxxx;
	}
}

void KS1dEigen::jacobian(const state_type & u,
    const scalar_type /*t*/,
	jacobian_type & jac)const{

  //evaluate jacobian
  if (jac.rows() == 0 || jac.cols()==0 ){
    jac.resize(u.size(), u.size());
  }
  tripletList.clear();
  tripletList.push_back( Tr( 0, 0, -dxInv_*u(0)) );

  scalar_type dxInv2 = dxInv_ * dxInv_;
  scalar_type dxInv4 = dxInv2 * dxInv2;

  for (ui_t i=0; i<Nnode_; ++i){

	  scalar_type Jmm = 0.0;
	  scalar_type Jm = 0.0;
	  scalar_type Jc = 0.0;
	  scalar_type Jp = 0.0;
	  scalar_type Jpp = 0.0;

		// Boundary conditions
	  scalar_type up = (i < Nnode_ - 1) ? u(i+1) : 0;//u[0];
	  scalar_type um = (i > 0) ? u(i-1) : 0; //u[N_GRID-1];

	  //du/dx
	  // -c * 0.5 * dxInv_ * (up - um);
	  Jm += mu_(0) * 0.5 * dxInv_ * um;
	  Jp += -1 * mu_(0) * 0.5 * dxInv_ * up;

      //du2/dx
	  // -0.5 * dxInv_ * (up*vp - um*vm) ;
	  Jm += 0.5 * dxInv_ * um;
	  Jp += -0.5 * dxInv_ * up;

	  //d2u/dx2
	  // -1 * dxInv_ * dxInv_ * (up + um - 2 * u(i)) ;
	  Jm += -1 * dxInv2;
	  Jc += 2 * dxInv2;
	  Jp += -1 * dxInv2;

	  //d4u/dx4
	  // -1 * dxInv_^4 * (upp - 4*up + 6 * u(i) -4* um + umm) ;
	  if (i > 1) Jmm += - mu_(1) *dxInv4;
	  Jm += 4.0 *  mu_(1) * dxInv4;

	  // If statements to set boundary conditions
	  if ((i<2) || ( i>Nnode_-3))
	  {
		  Jc += -7 *  mu_(1) * dxInv4;
	  }
	  else
	  {
		  Jc += -6 *  mu_(1) * dxInv4;
	  }

	  Jp += 4.0 * dxInv4;
	  if (i < Nnode_ - 2) Jpp += - mu_(1) *dxInv4;

	  // Store in triplet list
	  if (i > 1) tripletList.push_back( Tr( i, i-2, Jmm) );
	  if (i > 0) tripletList.push_back( Tr( i, i-1, Jm) );
      tripletList.push_back( Tr( i, i, Jc ) );
	  if (i < Nnode_-1) tripletList.push_back( Tr( i, i+1, Jp) );
	  if (i < Nnode_-2) tripletList.push_back( Tr( i, i+2, Jpp) );

  }
  jac.setFromTriplets(tripletList.begin(), tripletList.end());
}

}} //namespace pressio::apps
