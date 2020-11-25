
#ifndef PRESSIO_TESTS_SOLVERS_FROZEN_JAC_SHARED_TYPES_HPP_
#define PRESSIO_TESTS_SOLVERS_FROZEN_JAC_SHARED_TYPES_HPP_

#include "pressio_containers.hpp"

using eig_mat = Eigen::Matrix<double, -1, -1>;
using eig_vec = Eigen::VectorXd;
using vec_type = pressio::containers::Vector<eig_vec>;
using mat_type = pressio::containers::DenseMatrix<eig_mat>;

#endif
