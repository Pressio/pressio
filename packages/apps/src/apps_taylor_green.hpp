
#ifndef APPS_TAYLOR_GREEN_HPP_
#define APPS_TAYLOR_GREEN_HPP_

#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cassert>

#include "apps_ConfigDefs.hpp"
#include "apps_constants.hpp"
#include "apps_mpi_square_grid.hpp"
#include "apps_square_periodic_cc_mesh.hpp"
#include "apps_cc_vort_stream_square_periodic_second_order.hpp"

namespace apps{


class TaylorGreen2d
  : public apps::impl::VortStCCSqPeriodicSecOrd{

public:
  using base_t = apps::impl::VortStCCSqPeriodicSecOrd;
  using scalar_type = double;
  using state_type = Epetra_Vector;
  using space_residual_type = Epetra_Vector;

  struct TGVort{
    double operator()(double x, double y, double t){
      return -2.0*std::cos(x)*std::cos(y)*std::exp(-2.*t/Re_);
    }
  };
  struct TGStream{
    double operator()(double x, double y, double t){
      return -std::cos(x)*std::cos(y)*std::exp(-2.*t/Re_);}
  };
  
public:  
  TaylorGreen2d(int Ncells, Epetra_MpiComm * comm)
    : base_t(comm, Re_), N_(Ncells){}
  ~TaylorGreen2d() = default; 

public:
  void setup(){
    meshInfo_.setup(xL_, xR_, N_);
    base_t::setup();
  }// end
  //-----------------------------------
  
  void fixRHSForPoissonSolve(double time) final{
    if (mpiInfo_.myR_==0) {
      (*b_)[0] = TGStream()((*meshInfo_.xx_)[0], (*meshInfo_.yy_)[0], time);
    }
  }//
  //-----------------------------------

  void initializeVorticity() final{
    for (int i=0; i<meshInfo_.lNsq_; i++)
      (*omeNoGh_)[i]=TGVort()((*meshInfo_.xx_)[i],(*meshInfo_.yy_)[i], 0.0);
  }//
  //-----------------------------------
  
  double computeError(const Epetra_Vector & myF, double time, std::string f){
    double err = 0.0, tmp = 0.0;
    for (int i=0; i<meshInfo_.lNsq_; i++){
      if (f == "sf")
  	tmp = std::abs(myF[i] - TGStream()((*meshInfo_.xx_)[i], (*meshInfo_.yy_)[i], time));
      else if (f == "vo")
  	tmp = std::abs(myF[i] - TGVort()((*meshInfo_.xx_)[i], (*meshInfo_.yy_)[i], time));

      err += tmp*tmp;
      // if (myR_==0){
      // 	std::cout << std::setprecision(14)
      // 		  <<  (*xx_)[i] << " "
      // 		  <<  (*yy_)[i] << " "
      // 		  << (*UU_)[i] << " "
      // 		  << vel_U((*xx_)[i],(*yy_)[i],time)
      // 		  << std::endl;
      // }
    }
    double gErr = 0.0;
    mpiInfo_.comm_->SumAll(&err, &gErr, mpiInfo_.nProc_);
    return gErr;
  }//end 


  // void run(double dt, int Nsteps = 1)
  // {
  //   double t = 0.0;
  //   initializeVorticity();
  //   // time loop
  //   for (int iStep=1; iStep<=Nsteps; iStep++){      
  //     residual( *omeNoGhK1_, *omeNoGh_, t);
  //     double errPSI = computeError(*psiNoGh_, t, "sf");
  //     if (mpiInfo_.myR_==0){
  //     	std::cout << " step: " << iStep
  //     		  << " errSF = " << std::setprecision(6)
  //     		  << std::sqrt(errPSI/(meshInfo_.N_ * meshInfo_.N_))
  //     		  << std::endl;
  //     }
  //     for (int i=0; i<meshInfo_.lNsq_; i++)
  //     	(*omeNoGh_)[i] += dt * (*omeNoGhK1_)[i];
  //     t = static_cast<double>(iStep) * dt;
  //     double errVOR = computeError(*omeNoGh_, t, "vo");
  //     if (mpiInfo_.myR_==0){
  //     	std::cout << " step: " << iStep
  //     		  << " errVO = " << std::setprecision(6)
  //     		  << std::sqrt(errVOR/(meshInfo_.N_*meshInfo_.N_))
  //     		  << std::endl;
  //     }
      
  //     // *yn_ = *omeNoGh_;      
  //     // //RKs1
  //     // residual( *omeNoGhK1_, *yn_, t);
  //     // for (int i=0; i<lNsq_; i++)
  //     // 	(*omeNoGh_)[i] = (*yn_)[i] + (2./3.) * dt * (*omeNoGhK1_)[i];
  //     // t += (2./3.) * dt;
  //     // //RKs2
  //     // residual( *omeNoGhK2_, *omeNoGh_, t);
  //     // for (int i=0; i<lNsq_; i++)
  //     // 	(*omeNoGh_)[i] = (*yn_)[i]
  //     //  + dt * ( 0.25*(*omeNoGhK1_)[i] + (3./4.)*(*omeNoGhK2_)[i] );
  //     // t = static_cast<double>(iStep) * dt;
      
  //     // if (iStep % 10 == 0)
  //     // 	  printField("C_" + std::to_string(iStep) + ".txt", *omeNoGh_);
  //     // if (myR_==0){
  //     // 	std::cout << " done with step: " << iStep
  //     // 		  << " t = " << t << std::endl;
  //     // }      
  //     // double errVOR = computeError(*omeNoGh_, t, "vo");
  //     // if (myR_==0){
  //     // 	std::cout << " step: " << iStep
  //     // 		  << " errVO = " << std::setprecision(6)
  //     // 		  << std::sqrt(errVOR/(N_*N_))
  //     // 		  << std::endl;
  //     // }
  //   }

  // }//end 
  // //-----------------------------------


private:
  int N_; // # of cells along each axis
  const scalar_type xL_ = 0.0; //left side of domain 
  const scalar_type xR_ = 2.*apps::impl::PI_; // right side of domain
  static constexpr double Re_ = 10000.;  
  // apps::impl::MpiSquareGrid mpiInfo_;
  // apps::impl::SquarePeriodicCCMesh grid_;

}; //end class
  
} //end namespace apps
#endif 
