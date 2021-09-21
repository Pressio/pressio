/*
//@HEADER
// ************************************************************************
//
// tutorial1.cc
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <array>
#include <Eigen/Core>

#include "pressio/type_traits.hpp"
#include "pressio/ops.hpp"

struct MyCustomVector
{
  MyCustomVector(std::size_t ext) : d_(ext){}

  double & operator()(int i){ return d_[i]; }
  const double & operator()(int i)const { return d_[i]; }

  std::size_t extent(int k)const { return (k==0) ? d_.size() : 0; }

  void fill(double value){
    std::for_each(d_.begin(), d_.end(), [](double & v){ v= 0.; });
  }

private:
  std::vector<double> d_ = {};
};

struct MyCustomMatrix
{
  MyCustomMatrix(std::size_t nr, std::size_t nc)
    : num_rows_(nr), num_cols_(nc), d_(nr*nc){}

  std::size_t extent(int k)const { return (k==0) ? num_rows_ : num_cols_; }

  double & operator()(int i, int j){ return d_[num_cols_*i+j]; }
  const double & operator()(int i, int j) const { return d_[num_cols_*i+j]; }

  void fill(double value){
    std::for_each(d_.begin(), d_.end(), [=](double & v){ v=value; });
  }

private:
  std::size_t num_rows_ = {};
  std::size_t num_cols_ = {};
  std::vector<double> d_ = {};
};

struct MyRosenbrockSystem
{
  using scalar_type   = double;
  using state_type    = Eigen::VectorXd;
  using residual_type = MyCustomVector;
  using jacobian_type = MyCustomMatrix;

  residual_type createResidual() const{ return residual_type(6);   }
  jacobian_type createJacobian() const{ return jacobian_type(6, 4);}

  void residual(const state_type& x, residual_type & res) const
  {
    const auto & x1 = x(0);
    const auto & x2 = x(1);
    const auto & x3 = x(2);
    const auto & x4 = x(3);
    res(0) = 10.*(x4 - x3*x3);
    res(1) = 10.*(x3 - x2*x2);
    res(2) = 10.*(x2 - x1*x1);
    res(3) = (1.-x1);
    res(4) = (1.-x2);
    res(5) = (1.-x3);
  }

  void jacobian(const state_type & x, jacobian_type & JJ) const
  {
    const auto & x1 = x(0);
    const auto & x2 = x(1);
    const auto & x3 = x(2);
    JJ.fill(0.);

    JJ(0,2) = -20.*x3;
    JJ(0,3) = 10.;
    JJ(1,1) = -20.*x2;
    JJ(1,2) = 10.;
    JJ(2,0) = -20.*x1;
    JJ(2,1) = 10.;
    JJ(3,0) = -1.;
    JJ(4,1) = -1.;
    JJ(5,2) = -1.;
  }
};

namespace pressio{
template<> struct Traits<MyCustomVector>{
  using scalar_type = double;
};

template<> struct Traits<MyCustomMatrix>{
  using scalar_type = double;
};

namespace ops{
MyCustomVector clone(const MyCustomVector & src){ return src; }
MyCustomMatrix clone(const MyCustomMatrix & src){ return src; }

void set_zero(MyCustomVector & o){ o.fill(0); }
void set_zero(MyCustomMatrix & o){ o.fill(0); }

double norm2(const MyCustomVector & v){
  double norm{0};
  for (std::size_t i=0; i<v.extent(0); ++i){
    norm += v(i)*v(i);
  }
  return std::sqrt(norm);
}

template<class HessianType>
void product(pressio::transpose, pressio::nontranspose,
	     const double alpha,
	     const MyCustomMatrix & A,
	     const double beta,
	     HessianType & H)
{
  for (std::size_t i=0; i<A.extent(1); ++i){
    for (std::size_t j=0; j<A.extent(1); ++j)
    {
      H(i,j) *= beta;
      for (std::size_t k=0; k<A.extent(0); ++k){
	H(i,j) += alpha * A(k,i) * A(k,j);
      }
    }
  }
}

template<class GradientType>
void product(pressio::transpose,
	     const double alpha,
	     const MyCustomMatrix & A,
	     const MyCustomVector & b,
	     const double beta,
	     GradientType & g)
{
  for (int i=0; i<g.rows(); ++i){
    g(i) *= beta;
    for (std::size_t k=0; k<A.extent(0); ++k){
      g(i) += alpha * A(k,i) * b(k);
    }
  }
}

template<class HessianType>
HessianType product(pressio::transpose, pressio::nontranspose,
	     double alpha,
	     const MyCustomMatrix & A)
{
  HessianType H(A.extent(1), A.extent(1));
  product(pressio::transpose(), pressio::nontranspose(), alpha, A, 0, H);
  return H;
}

void update(MyCustomVector & v, double a, const MyCustomVector & v1, double b)
{
  for (std::size_t i=0; i<v.extent(0); ++i){
    v(i) = v(i)*a + b*v1(i);
  }
}

void scale(MyCustomVector & v, double factor){
  for (std::size_t i=0; i<v.extent(0); ++i){
    v(i) = v(i)*factor;
  }
}
}}//end namespace pressio::ops


#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear.hpp"

int main()
{
  namespace plog   = pressio::log;
  namespace pls    = pressio::linearsolvers;
  namespace pnonls = pressio::nonlinearsolvers;
  plog::initialize(pressio::logto::terminal);
  plog::setVerbosity({plog::level::info});

  using problem_t = MyRosenbrockSystem;
  problem_t problem;

  using state_t   = Eigen::VectorXd;
  state_t x(4);
  x[0] = -0.05; x[1] = 1.1; x[2] = 1.2; x[3] = 1.5;

  using hessian_t    = Eigen::MatrixXd;
  using lin_tag      = pls::direct::HouseholderQR;
  using lin_solver_t = pls::Solver<lin_tag, hessian_t>;
  lin_solver_t linSolver;

  auto gnSolver = pnonls::create_gauss_newton(problem, x, linSolver);
  gnSolver.setTolerance(1e-5);
  gnSolver.solve(problem, x);
  std::cout << std::setprecision(14) << x << std::endl;
  // check solution
  std::cout << "Computed solution: \n "
            << "[" << x(0) << " " << x(1) << " " << x(2) << " " << x(3) << " " << "] \n"
            << "Expected solution: \n "
            << "[1.0000000156741, 0.99999999912477, 0.99999999651993, 0.99999998889888]"
            << std::endl;


  std::vector<double> gold = {
    1.00000001567414e+00,
    9.99999999124769e-01,
    9.99999996519930e-01,
    9.99999988898883e-01};

  std::string sentinel = "PASSED";
  const auto e1 = std::abs(x(0) - gold[0]);
  const auto e2 = std::abs(x(1) - gold[1]);
  const auto e3 = std::abs(x(2) - gold[2]);
  const auto e4 = std::abs(x(3) - gold[3]);
  if (e1>1e-6 or e2>1e-6  or e3>1e-6 or e4>1e-6){
    sentinel = "FAILED";
  }
  std::cout << sentinel << std::endl;

  plog::finalize();
  return 0;
}
