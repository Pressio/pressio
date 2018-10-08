
#ifndef ODE_REF_APPS_FOR_TESTING_HPP_
#define ODE_REF_APPS_FOR_TESTING_HPP_

namespace rompp{ namespace ode{ namespace testing{ 

      
//************************************************
struct fakeAppForTraitsForExp{
  using scalar_type = double;
  using state_type = std::vector<double>;
  using residual_type = std::vector<double>;

  void residual(const state_type & y,
		residual_type & R,
		scalar_type t) const{
  };
  residual_type residual(const state_type & y, scalar_type t )const{
    residual_type R;
    return R;
  };
};
//************************************************
//************************************************

      
struct refAppEigen{

  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;

  void residual(const state_type & y,
		residual_type & R,
		scalar_type t) const{
    auto sz = y.size();
    for (decltype(sz) i=0; i<sz; i++)
      R[i] = y[i];
  };
  
  residual_type residual(const state_type & y,
			 scalar_type t) const{
    residual_type R(y);
    residual(y, R, t);
    return R;
  };

};
//************************************************
//************************************************

      
struct refAppForImpEigen{

  /*
    dy
    -- = -10*y
    dt 

    y(0) = 1; 
    y(1) = 2; 
    y(2) = 3;

    for a given time-step, dt, we have: 
       Euler backward yields: 
          y_n+1 = y_n/(1+10*dt)
   */
  
  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  state_type y0;

public:  
  refAppForImpEigen(){
    y0.resize(3);
    y0 << 1., 2., 3.;
  }

  void residual(const state_type & y, residual_type & R,
		scalar_type t) const{
    auto sz = y.size(); assert(sz==3);
    R = -10. * y;
  };
  //--------------------------------------------

  residual_type residual(const state_type & y,
			 scalar_type t) const{
    residual_type R(y);
    assert(R.size()==3);
    residual(y, R, t);
    return R;
  };
  //--------------------------------------------
  
  void jacobian(const state_type & y,
		jacobian_type & JJ,
		scalar_type t) const{
    assert( JJ.rows() == 3 ); assert( JJ.cols() == 3 );

    typedef Eigen::Triplet<scalar_type> Tr;
    std::vector<Tr> tripletList;
    tripletList.push_back( Tr( 0, 0, -10.) );
    tripletList.push_back( Tr( 1, 1, -10.) );
    tripletList.push_back( Tr( 2, 2, -10.) );
    JJ.setFromTriplets(tripletList.begin(), tripletList.end());
  };
  //--------------------------------------------
  
  jacobian_type jacobian(const state_type & y,
			 scalar_type t) const{
    jacobian_type JJ(3,3);
    jacobian(y, JJ, t);
    return JJ;
  };
  //--------------------------------------------

  void analyticAdvanceBackEulerNSteps(double dt, int n){
    double den = 1.0 + 10.*dt;
    for (int i=1; i!=n+1; i++){
      y0[0] = y0[0]/den;
      y0[1] = y0[1]/den;
      y0[2] = y0[2]/den;
    }
  };

};//end app refAppForImpEigen
//************************************************
//************************************************

}}} // namespace rompp::ode::testing
#endif
