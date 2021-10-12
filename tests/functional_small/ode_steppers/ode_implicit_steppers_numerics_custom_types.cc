
#include <gtest/gtest.h>

using ScalarType = double;
using VectorType = std::vector<ScalarType>;
using MatrixType = std::vector<std::vector<ScalarType>>;

namespace
{

void resize_jacobian(MatrixType & M){
  M[0].resize(3);
  M[1].resize(3);
  M[2].resize(3);
}

struct MyApp{
  using scalar_type   = ScalarType;
  using state_type    = VectorType;
  using velocity_type = VectorType;
  using jacobian_type = MatrixType;

public:
  void velocity(const state_type & y,
	   	          scalar_type t, 
                velocity_type & f) const
  {
    f[0] = 1.;
    f[1] = 1.;
    f[2] = 1.;
  };

  void jacobian(const state_type & yIn,
                scalar_type t, 
                jacobian_type & J) const
  {
    for (auto & it : J){
      for (auto & it2 : it){
        it2 = 1.;
      }
    }
  }

  velocity_type createVelocity() const{
    velocity_type R(3);
    return R;
  };

  jacobian_type createJacobian() const
  {
    jacobian_type J(3);
    resize_jacobian(J);
    return J;

  }
};
} // end namespace

#include "pressio/type_traits.hpp"

namespace pressio{ 

template<> struct Traits<VectorType>{
  using scalar_type = ScalarType;
};

template<> struct Traits<MatrixType>{
  using scalar_type = ScalarType;
};

namespace ops{
void deep_copy(VectorType & dest, const VectorType & from){
  dest = from;
}

VectorType clone(const VectorType & src){
  return VectorType(src.size());
}

void scale(MatrixType & M, ScalarType factor){
  for (auto & it : M){
    for (auto & it2 : it){
      it2 *= factor;
    }
  }
}

void add_to_diagonal(MatrixType & M, ScalarType value)
{
  for (int i=0; i<3; ++i){
    M[i][i] += value;
  }
}

void update(VectorType & v,        const ScalarType a,
		        const VectorType & v1, const ScalarType b)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v1[i];
  }
}

void update(VectorType & v,        const ScalarType a,
            const VectorType & v0, const ScalarType b,
            const VectorType & v1, const ScalarType c)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v0[i] + c*v1[i];
  }
}

void update(VectorType & v, const ScalarType a,
    const VectorType & v1, const ScalarType b,
    const VectorType & v2, const ScalarType c,
    const VectorType & v3, const ScalarType d)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v1[i] + c*v2[i] + d*v3[i];
  }
}

void update(VectorType & v, const ScalarType a,
    const VectorType & v1, const ScalarType b,
    const VectorType & v2, const ScalarType c,
    const VectorType & v3, const ScalarType d,
    const VectorType & v4, const ScalarType e)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
  }
}
}} //end namespace pressio::ops

#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"

struct MySolverForBDF1
{
  template<class StepperType>
  void solve(StepperType & stepperObj, VectorType & x)
  {
    VectorType R(3);
    MatrixType J(3);
    resize_jacobian(J);

    EXPECT_DOUBLE_EQ(x[0], 1.);
    EXPECT_DOUBLE_EQ(x[1], 2.);
    EXPECT_DOUBLE_EQ(x[2], 3.);

    stepperObj.residual(x, R);
    EXPECT_DOUBLE_EQ(R[0], -10.);
    EXPECT_DOUBLE_EQ(R[1], -10.);
    EXPECT_DOUBLE_EQ(R[2], -10.);

    stepperObj.jacobian(x, J);
    for (int i=0; i<3; ++i){
      for (int j=0; j<3; ++j){
        if (i==j) {
          EXPECT_DOUBLE_EQ(J[i][j], -9.);
        } 
        else{
          EXPECT_DOUBLE_EQ(J[i][j], -10.);
        }
      }
    }

    x[0]=1.1;
    x[1]=2.2;
    x[2]=3.3;
  }
};

TEST(ode, implicit_bdf1_custom_types)
{
  using namespace pressio;
  using app_t	   = MyApp;
  using state_t = typename app_t::state_type;
  app_t appObj;

  state_t y(3);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;
  auto stepperObj = ode::create_bdf1_stepper(y,appObj);

  ScalarType dt = 10.;
  MySolverForBDF1 solver;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 1, solver);
  EXPECT_DOUBLE_EQ( y[0], 1.1);
  EXPECT_DOUBLE_EQ( y[1], 2.2);
  EXPECT_DOUBLE_EQ( y[2], 3.3);
}


struct MySolverForBDF2
{
  int count_ = 0;

  template<class StepperType>
  void solve(StepperType & stepperObj, VectorType & x)
  {
    ++count_;

    VectorType R(3);
    MatrixType J(3);
    resize_jacobian(J);

    if (count_==1){
      EXPECT_DOUBLE_EQ(x[0], 1.);
      EXPECT_DOUBLE_EQ(x[1], 2.);
      EXPECT_DOUBLE_EQ(x[2], 3.);

      stepperObj.residual(x, R);
      EXPECT_DOUBLE_EQ(R[0], -10.);
      EXPECT_DOUBLE_EQ(R[1], -10.);
      EXPECT_DOUBLE_EQ(R[2], -10.);

      stepperObj.jacobian(x, J);
      for (int i=0; i<3; ++i){
        for (int j=0; j<3; ++j){
          if (i==j) {
            EXPECT_DOUBLE_EQ(J[i][j], -9.);
          } 
          else{
            EXPECT_DOUBLE_EQ(J[i][j], -10.);
          }
        }
      }
    }

    if (count_==2){
      stepperObj.residual(x, R);
      EXPECT_DOUBLE_EQ(R[0], 1.-(4/3.)*1. + (1./3)*1 - (2./3.)*10.);
      EXPECT_DOUBLE_EQ(R[1], 2.-(4/3.)*2. + (1./3)*2 - (2./3.)*10.);
      EXPECT_DOUBLE_EQ(R[2], 3.-(4/3.)*3. + (1./3)*3 - (2./3.)*10.);

      stepperObj.jacobian(x, J);
      for (int i=0; i<3; ++i){
        for (int j=0; j<3; ++j){
          if (i==j) {
            EXPECT_DOUBLE_EQ(J[i][j], -(2./3)*10.*1. + 1.);
          } 
          else{
            EXPECT_DOUBLE_EQ(J[i][j], -(2./3)*10.*1.);
          }
        }
      }
      x[0]=1.1;
      x[1]=2.2;
      x[2]=3.3;
    }
  }
};

TEST(ode, implicit_bdf2_custom_types)
{
  using namespace pressio;
  using app_t    = MyApp;
  using state_t = typename app_t::state_type;
  app_t appObj;

  state_t y(3);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;
  auto stepperObj = ode::create_bdf2_stepper(y,appObj);

  ScalarType dt = 10.;
  MySolverForBDF2 solver;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 2, solver);
  EXPECT_DOUBLE_EQ( y[0], 1.1);
  EXPECT_DOUBLE_EQ( y[1], 2.2);
  EXPECT_DOUBLE_EQ( y[2], 3.3);
}


struct MySolverForCrankNic
{
  int count_ = 0;

  template<class StepperType>
  void solve(StepperType & stepperObj, VectorType & x)
  {
    ++count_;

    VectorType R(3);
    MatrixType J(3);
    resize_jacobian(J);

    if (count_==1){
      EXPECT_DOUBLE_EQ(x[0], 1.);
      EXPECT_DOUBLE_EQ(x[1], 2.);
      EXPECT_DOUBLE_EQ(x[2], 3.);

      stepperObj.residual(x, R);
      EXPECT_DOUBLE_EQ(R[0], -0.5*10.*2.);
      EXPECT_DOUBLE_EQ(R[1], -0.5*10.*2.);
      EXPECT_DOUBLE_EQ(R[2], -0.5*10.*2.);

      stepperObj.jacobian(x, J);
      for (int i=0; i<3; ++i){
        for (int j=0; j<3; ++j){
          if (i==j) {
            EXPECT_DOUBLE_EQ(J[i][j], -4.);
          } 
          else{
            EXPECT_DOUBLE_EQ(J[i][j], -5.);
          }
        }
      }
    }

    x[0]=1.1;
    x[1]=2.2;
    x[2]=3.3;
  }
};

TEST(ode, implicit_crank_nicolson_custom_types)
{
  using namespace pressio;
  using app_t    = MyApp;
  using state_t = typename app_t::state_type;
  app_t appObj;

  state_t y(3);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;
  auto stepperObj = ode::create_cranknicolson_stepper(y,appObj);

  ScalarType dt = 10.;
  MySolverForCrankNic solver;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 1, solver);
  EXPECT_DOUBLE_EQ( y[0], 1.1);
  EXPECT_DOUBLE_EQ( y[1], 2.2);
  EXPECT_DOUBLE_EQ( y[2], 3.3);
}
