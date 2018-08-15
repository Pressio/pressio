
#include <gtest/gtest.h>

#include <iostream>

#include "experimental/solvers_policy_nonlinear_iterative.hpp"
#include "matrix/concrete/core_matrix_dense_serial_eigen.hpp"
#include "matrix/concrete/core_matrix_sparse_serial_eigen.hpp"
#include "vector/concrete/core_vector_serial_eigen.hpp"


struct ValidSystem {

    // Matrix typedefs
    using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_w_t = core::Matrix<matrix_n_t>;

    // Vector typedefs
    using vector_n_t = Eigen::VectorXd;
    using vector_w_t = core::Vector<vector_n_t>;

    typedef vector_w_t vector_type;
    typedef matrix_w_t matrix_type;


    void residual(const vector_w_t& x, vector_w_t& res) {
      res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
      res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
    }


    auto residual(const vector_w_t& x) {
      vector_w_t res(2);
      res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
      res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
      return res;
    }

    
    void jacobian(const vector_w_t& x, matrix_w_t& jac) {
      jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
      jac.data()->coeffRef(0, 1) =  1.0;
      jac.data()->coeffRef(1, 0) = -1.0;
      jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1]; 
    }


    auto jacobian(const vector_w_t& x) {
      matrix_w_t jac(2, 2);
      jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
      jac.data()->coeffRef(0, 1) =  1.0;
      jac.data()->coeffRef(1, 0) = -1.0;
      jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1]; 
      return jac;
    }
};


TEST(solvers_nonlinear_iterative_newtonraphson, solvers_nonlinear_iterative_newtonraphsonPolicyCompatibleSystemVector)
{
  using namespace solvers;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = core::Vector<vector_n_t>;

  typedef vector_w_t vector_type;

  vector_w_t b;

  auto val = SolversNonLinearIterativeNewtonRaphsonPolicy::solve(ValidSystem{}, b);

  EXPECT_EQ(val, 0);
}

TEST(solvers_nonlinear_iterative_newtonraphson, solvers_nonlinear_iterative_newtonraphsonPolicyIncompatibleSystemVector)
{
  using namespace solvers;

  int b;

  ASSERT_DEATH(SolversNonLinearIterativeNewtonRaphsonPolicy::solve(ValidSystem{}, b), "Error: the type of the RHS vector is not compatible with the provided nonlinear system");
}