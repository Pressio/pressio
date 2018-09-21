
#ifndef APPS_SHEAR_LAYER_HPP_
#define APPS_SHEAR_LAYER_HPP_

#include <cmath>
#include "apps_ConfigDefs.hpp"
#include "apps_constants.hpp"
#include "apps_mpi_square_grid.hpp"
#include "apps_square_periodic_cc_mesh.hpp"
#include "apps_cc_vort_stream_square_periodic_second_order.hpp"

namespace apps{


class ShearLayer
  : public apps::impl::VortStCCSqPeriodicSecOrd{

  static constexpr double PI = apps::impl::PI_;
  static constexpr double PIh = apps::impl::PI_*0.5;
  static constexpr double PI3h = apps::impl::PI_*1.5;

public:
  using base_t = apps::impl::VortStCCSqPeriodicSecOrd;
  using scalar_type = double;
  using state_type = Epetra_Vector;
  using space_residual_type = Epetra_Vector;

  struct SHLAVort{
    double operator()(double x, double y, double t){
      double vortX = delta_*std::cos(x);
      if (y <= PI){
	double y2 = sig_*(y - PIh);
	double den = std::cosh(y2);
	return vortX - sig_/(den*den);
      }
      else {
	double y2 = sig_*(PI3h-y);
	double den = std::cosh(y2);
	return vortX + sig_/(den*den);
      }
    }
  };

  struct SHLAUU{
    double operator()(double x, double y, double t)
    {
      if (y <= PI){
	double y2 = sig_*(y - PIh);
	return std::tanh(y2);
      }
      else {
	double y2 = sig_*(PI3h-y);
	return std::tanh(y2);
      }
    }
  };

  struct SHLAVV{
    double operator()(double x, double y, double t){
      return delta_*std::sin(x);
    }
  };
  
public:  
  ShearLayer(int Ncells, Epetra_MpiComm * comm)
    : base_t(comm, Re_), N_(Ncells){}
  ~ShearLayer() = default; 

public:
  void setup(){
    meshInfo_.setup(xL_, xR_, N_);
    base_t::setup();
  }// end
  //-----------------------------------
  
  void fixRHSForPoissonSolve(double time) final{
    if (mpiInfo_.myR_==0) {
      (*b_)[0] = 0.;//TGStream()((*meshInfo_.xx_)[0], (*meshInfo_.yy_)[0], time);
    }
  }//
  //-----------------------------------

  void initializeVorticity() final{
    for (int i=0; i<meshInfo_.lNsq_; i++)
      (*omeNoGh_)[i]=SHLAVort()((*meshInfo_.xx_)[i],(*meshInfo_.yy_)[i], 0.0);
  }//
  //-----------------------------------

  double computeError(const Epetra_Vector & myF, double time, std::string f){
    double err = 0.0, tmp = 0.0;
    for (int i=0; i<meshInfo_.lNsq_; i++){
      double x = (*meshInfo_.xx_)[i];
      double y = (*meshInfo_.yy_)[i];
      
      if (f == "uu")
  	tmp = std::abs(myF[i] -SHLAUU()(x,y,time));
      else if (f == "vv")
  	tmp = std::abs(myF[i] -SHLAVV()(x,y,time));

      err += tmp*tmp;
      // if (mpiInfo_.myR_==0){
      // 	std::cout << std::setprecision(14)
      // 		  << x << " "
      // 		  << y << " "
      // 		  << myF[i] << " "
      // 		  << TGUU()(x,y,time)
      // 		  << std::endl;
      // }
    }
    double gErr = 0.0;
    mpiInfo_.comm_->SumAll(&err, &gErr, mpiInfo_.nProc_);
    return gErr;
  }//end 
  //------------------------------------------------

  
  void run(double dt, int Nsteps = 1)
  {
    double t = 0.0;
    initializeVorticity();

    // time loop
    for (int iStep=1; iStep<=Nsteps; iStep++)
    { 
      *yn_ = *omeNoGh_;      

      //RKs1
      residual( *omeNoGhK1_, *yn_, t);
      for (int i=0; i<meshInfo_.lNsq_; i++)
      	(*omeNoGh_)[i] = (*yn_)[i] + (2./3.) * dt * (*omeNoGhK1_)[i];
      t += (2./3.) * dt;

      //RKs2
      residual( *omeNoGhK2_, *omeNoGh_, t);
      for (int i=0; i<meshInfo_.lNsq_; i++)
      	(*omeNoGh_)[i] = (*yn_)[i]
       + dt * ( 0.25*(*omeNoGhK1_)[i] + (3./4.)*(*omeNoGhK2_)[i] );
      t = static_cast<double>(iStep) * dt;
      
      // // residual( *omeNoGhK1_, *omeNoGh_, t);
      // // for (int i=0; i<meshInfo_.lNsq_; i++)
      // // 	(*omeNoGh_)[i] += dt * (*omeNoGhK1_)[i];
      // // t = static_cast<double>(iStep) * dt;
      // double err = computeError(*UU_, t, "uu");
      // if (mpiInfo_.myR_==0){
      // 	std::cout << " step: " << iStep
      // 		  << " err = " << std::setprecision(6)
      // 		  << std::sqrt(err/(N_*N_))
      // 		  << std::endl;
      // }      
      if (iStep % 100 == 0)
      	printField("C_" + std::to_string(iStep)
      		   + "_P" + std::to_string(mpiInfo_.myR_) + ".txt", *omeNoGh_);
      if (mpiInfo_.myR_==0){
      	std::cout << " done with step: " << iStep
      		  << " t = " << t << std::endl;
      }      
    }

  }//end 
  //-----------------------------------
  
private:

  // see for reference:  http://dx.doi.org/10.1016/j.jcp.2012.09.005
    
  int N_; // # of cells along each axis
  const scalar_type xL_ = 0.0; //left side of domain 
  const scalar_type xR_ = 2.*apps::impl::PI_; // side of domain
  static constexpr double Re_ = 10000.; 
  static constexpr double sig_ = 15./apps::impl::PI_;
  static constexpr double delta_ = 0.05;

}; //end class
  
} //end namespace apps
#endif 
