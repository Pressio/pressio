
#ifndef ODE_REF_APPS_FOR_TESTING_HPP_
#define ODE_REF_APPS_FOR_TESTING_HPP_

namespace rompp{ namespace ode{ namespace testing{

//************************************************
struct fakeAppForTraitsForExp{
  using scalar_type = double;
  using state_type = std::vector<double>;
  using velocity_type = std::vector<double>;

  void velocity(const state_type & y,
		velocity_type & R,
		scalar_type t) const{
  };
  velocity_type velocity(const state_type & y, scalar_type t )const{
    velocity_type R;
    return R;
  };
};
//************************************************
//************************************************

struct refAppEigen{
  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using velocity_type = state_type;

  void velocity(const state_type & y,
		velocity_type & R,
		scalar_type t) const{
    auto sz = y.size();
    for (decltype(sz) i=0; i<sz; i++)
      R[i] = y[i];
  };

  velocity_type velocity(const state_type & y,
			 scalar_type t) const{
    velocity_type R(y);
    velocity(y, R, t);
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

       BDF2 yields:
          num = (4/3) * y_n - (1/3) * y_n-2
	  den = 1 + (20/3) * dt
          y_n+1 = num/den
   */
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  state_type y;
  state_type y_nm1;

public:
  refAppForImpEigen(){
    y.resize(3);
    y << 1., 2., 3.;
  }

  void velocity(const state_type & y, velocity_type & R,
		scalar_type t) const{
    assert(y.size()==3);
    R = -10. * y;
  };
  //--------------------------------------------

  velocity_type velocity(const state_type & y,
			 scalar_type t) const{
    velocity_type R(y);
    assert(R.size()==3);
    velocity(y, R, t);
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
      y_nm1 = y;
      y[0] = y[0]/den;
      y[1] = y[1]/den;
      y[2] = y[2]/den;
    }
  };
  //--------------------------------------------

  void analyticAdvanceBDF2NSteps(double dt, int n){
    double den = 1.0 + (20./3.)*dt;
    for (int i=1; i!=n+1; i++){
      double num1 = (4./3.) * y[0] - (1./3.) * y_nm1[0];
      double num2 = (4./3.) * y[1] - (1./3.) * y_nm1[1];
      double num3 = (4./3.) * y[2] - (1./3.) * y_nm1[2];

      y_nm1 = y;
      y[0] = num1/den;
      y[1] = num2/den;
      y[2] = num3/den;
    }
  };

  void analyticAdvanceRK4(double dt)
  {
    assert(dt==0.1);
    velocity_type k1(3), k2(3), k3(3), k4(3);
    //I did the math...
    k1 << -1., -2, -3.;
    k2 << -0.5, -1, -1.5;
    k3 << -0.75, -1.5, -2.25;
    k4 << -0.25, -0.5, -0.75;
    y += (1./6.) * (k1 + 2*k2 + 2*k3 + k4);
  };
  //--------------------------------------------

};//end app refAppForImpEigen
//************************************************
//************************************************

}}} // namespace rompp::ode::testing
#endif
