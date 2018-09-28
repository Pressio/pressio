#include <iostream>
#include "CORE_MATRIX"
#include "SOLVERS_NONLINEAR"

struct NonLinearLeastSquareSystem {

    // Matrix typedefs
    using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_w_t = rompp::core::Matrix<matrix_n_t>;

    // Vector typedefs
    using vector_n_t = Eigen::VectorXd;
    using vector_w_t = rompp::core::Vector<vector_n_t>;

    typedef vector_w_t vector_type;
    typedef matrix_w_t matrix_type;


    void residual(const vector_w_t& x, vector_w_t& res) const {
      for (int i = 0; i < 5; i++) {
        res[i] = x[0] * exp(x[1] * coeff[i]) - vansw[i];
      }
    }


    auto residual(const vector_w_t& x) const {
      vector_w_t res(5);
      this->residual(x, res);
      return res;
    }


    void jacobian(const vector_w_t& x, matrix_w_t& jac) const {
      for (int i = 0; i < 5; i++) {
        double expval = exp(x[1] * coeff[i]);
        jac.data()->coeffRef(i, 0) = expval;
        jac.data()->coeffRef(i, 1) = x[0]*coeff[i]*expval;
      }
    }


    auto jacobian(const vector_w_t& x) const {
      matrix_w_t jac(5, 2);
      this->jacobian(x, jac);
      return jac;
    }

    const double coeff[5] = {1.0, 2.0, 4.0, 5.0, 8.0};
    const double vansw[5] = {3.2939, 4.2699, 7.1749, 9.3008, 20.259};
};

int main() {

  // Namespaces
  using namespace rompp;
  using namespace rompp::solvers;

  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = core::Vector<vector_n_t>;

  auto solver = NonLinearSolvers::createNonLinearIterativeLeastSquareSolver<nonlinearleastsquare::LevenbergMarquardt, linear::Bicgstab>();

  vector_w_t x0(2);
  x0[0] = 2.50;
  x0[1] = 0.25;

  NonLinearLeastSquareSystem sys;
  auto x = solver.solve(sys, x0);

  std::cout << "The solution of the nonlinear system is: " << std::endl << *x.data() << std::endl;

  return 0;
}
