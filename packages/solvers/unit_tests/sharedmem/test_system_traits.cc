
#include <gtest/gtest.h>

#include <iostream>

#include "experimental/system_traits.hpp"
#include "CORE_ALL"
// #include "matrix/concrete/core_matrix_sparse_serial_eigen.hpp"
// #include "vector/concrete/core_vector_serial_eigen.hpp"


struct ValidSystemA {

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


struct ValidSystemB {

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

    
    void jacobian(const vector_w_t& x, matrix_w_t& jac) {
      jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
      jac.data()->coeffRef(0, 1) =  1.0;
      jac.data()->coeffRef(1, 0) = -1.0;
      jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1]; 
    }
};


struct InvalidSystemA {

    // Matrix typedefs
    using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_w_t = core::Matrix<matrix_n_t>;

    // Vector typedefs
    using vector_n_t = Eigen::VectorXd;
    using vector_w_t = core::Vector<vector_n_t>;

    typedef vector_w_t vector_type;
    typedef matrix_w_t matrix_type;
    
    void jacobian(const vector_w_t& x, matrix_w_t& jac) {
      jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
      jac.data()->coeffRef(0, 1) =  1.0;
      jac.data()->coeffRef(1, 0) = -1.0;
      jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1]; 
    }
};


TEST(system_traits, systemTraitsValidSystemATest)
{
  using namespace solvers::details;

  typedef system_traits<ValidSystemA> system_t;

  auto is_system = system_t::is_system;

  auto has_residual_callable_with_one_arg = system_t::has_residual_callable_with_one_arg;
  auto has_residual_callable_with_two_args = system_t::has_residual_callable_with_two_args;

  auto has_jacobian_callable_with_one_arg = system_t::has_jacobian_callable_with_one_arg;
  auto has_jacobian_callable_with_two_args = system_t::has_jacobian_callable_with_two_args;

  EXPECT_EQ(is_system, true);
  EXPECT_EQ(has_residual_callable_with_one_arg, true);
  EXPECT_EQ(has_residual_callable_with_two_args, true);
  EXPECT_EQ(has_jacobian_callable_with_one_arg, true);
  EXPECT_EQ(has_jacobian_callable_with_two_args, true);
}


TEST(system_traits, systemTraitsValidSystemBTest) 
{
  using namespace solvers::details;

  typedef system_traits<ValidSystemB> system_t;

  auto is_system = system_t::is_system;

  auto has_residual_callable_with_one_arg = system_t::has_residual_callable_with_one_arg;
  auto has_residual_callable_with_two_args = system_t::has_residual_callable_with_two_args;

  auto has_jacobian_callable_with_one_arg = system_t::has_jacobian_callable_with_one_arg;
  auto has_jacobian_callable_with_two_args = system_t::has_jacobian_callable_with_two_args;

  EXPECT_EQ(is_system, true);
  EXPECT_EQ(has_residual_callable_with_one_arg, false);
  EXPECT_EQ(has_residual_callable_with_two_args, true);
  EXPECT_EQ(has_jacobian_callable_with_one_arg, false);
  EXPECT_EQ(has_jacobian_callable_with_two_args, true);
}


TEST(system_traits, systemTraitsInvalidSystemATest) 
{
  using namespace solvers::details;

  typedef system_traits<InvalidSystemA> system_t;

  auto is_system = system_t::is_system;

  auto has_residual_callable_with_one_arg = system_t::has_residual_callable_with_one_arg;
  auto has_residual_callable_with_two_args = system_t::has_residual_callable_with_two_args;

  auto has_jacobian_callable_with_one_arg = system_t::has_jacobian_callable_with_one_arg;
  auto has_jacobian_callable_with_two_args = system_t::has_jacobian_callable_with_two_args;

  EXPECT_EQ(is_system, false);
  EXPECT_EQ(has_residual_callable_with_one_arg, false);
  EXPECT_EQ(has_residual_callable_with_two_args, false);
  EXPECT_EQ(has_jacobian_callable_with_one_arg, false);
  EXPECT_EQ(has_jacobian_callable_with_two_args, true);
}