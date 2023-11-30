
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/rom_linear.hpp"

template<class T>
using column_expr_t =
  decltype( pressio::column( std::declval<T &>(), 0) );

class MyFomMatrixOrLinOperator{
  using _basis_type = Eigen::MatrixXd;
  using _vec_type = Eigen::VectorXd;

  _vec_type vecToAdd_;

public:
  explicit MyFomMatrixOrLinOperator(bool justMatrix) : vecToAdd_(10)
  {
    if (justMatrix){
      vecToAdd_.setZero();
    }
    else{
      vecToAdd_.setConstant(22.);
    }
  }

  _basis_type createResultOfActionOn(const _basis_type & operand) const{
    return _basis_type(10, operand.cols());
  }

  _vec_type createResultOfActionOn(const _vec_type & operand) const{
    return _vec_type(10);
  }

  void apply(column_expr_t<const _basis_type> operand,
	     column_expr_t<_basis_type> result) const
  {
    applyImpl(operand.native(), result.native());
  }

  void apply(const _vec_type & operand,
	     _vec_type & result) const
  {
    applyImpl(operand, result);
  }

private:
  template<class T1, class T2>
  void applyImpl(const T1 & operand, T2 & result) const  {
    result = operand;
    result *= 2;
    result += vecToAdd_;
  }
};

class ConstFomVector{
public:
  using vector_type = Eigen::VectorXd;

  vector_type evaluate() const{
    vector_type v(10);
    v.setConstant(55.);
    return v;
  }
};

TEST(rom, linear_rom_matrix_operator)
{
  using basis_t = Eigen::MatrixXd;
  using vec_t = Eigen::VectorXd;

  basis_t lb = basis_t::Random(10,6);
  basis_t rb = basis_t::Random(10,3);
  vec_t shift(10);
  shift.setConstant(11.);

  MyFomMatrixOrLinOperator fomO(true);

  namespace promlin = pressio::rom::linear;

  // - Asmall     = lb^T A rb
  // - Axrefsmall = lb^T A xref
  auto result = promlin::create_reduced_matrix_operator(fomO, lb, rb, shift);
  promlin::export_ascii(result,
			"linear_rom_red_mat_A.txt",
			"linear_rom_red_shift_A.txt");

  auto [A,s,_] = result();

  // the action of the fom A on a vector is just to scale it by 2
  // so we know what the gold should be
  auto gold = lb.transpose() * 2. * rb;
  ASSERT_TRUE(gold.isApprox(A));

  auto goldShift = lb.transpose() * 2. * shift;
  ASSERT_TRUE(goldShift.isApprox(s));
}

TEST(rom, linear_rom_vector_operator)
{
  using basis_t = Eigen::MatrixXd;

  basis_t lb = basis_t::Random(10,6);

  ConstFomVector fomO;
  namespace promlin = pressio::rom::linear;
  auto result = promlin::create_reduced_vector_operator(fomO, lb);
  promlin::export_ascii(result, "linear_rom_red_vec.txt");

  const auto & v = result();

  Eigen::VectorXd gold(10);
  gold.setConstant(1);
  gold = lb.transpose() * 55. * gold;
  ASSERT_TRUE(gold.isApprox(v));
}

TEST(rom, linear_rom_linear_operator)
{
  using basis_t = Eigen::MatrixXd;
  using vec_t = Eigen::VectorXd;

  basis_t lb = basis_t::Random(10,6);
  basis_t rb = basis_t::Random(10,3);
  vec_t shift(10);
  shift.setConstant(11.);

  MyFomMatrixOrLinOperator fomO(false);
  namespace promlin = pressio::rom::linear;
  auto result = promlin::create_reduced_linear_operator(fomO, lb, rb, shift);
  promlin::export_ascii(result,
			"linear_rom_red_mat_B.txt",
			"linear_rom_red_shift_B.txt",
			"linear_rom_red_vec_B.txt");

  auto [A,s,f] = result();

  std::cout << f << '\n';
  std::cout << "-----------------\n";
  Eigen::VectorXd goldForc(10);
  goldForc.setConstant(1);
  goldForc = lb.transpose() * 22. * goldForc;
  std::cout << goldForc << '\n';
  ASSERT_TRUE(goldForc.isApprox(f));

  // if everything is done correctly, then these should hold:
  auto gold = lb.transpose() * 2. * rb;
  ASSERT_TRUE(gold.isApprox(A));

  Eigen::VectorXd goldShift(10);
  goldShift.setConstant(0.);
  Eigen::VectorXd v2(10);
  v2.setConstant(44.);

  goldShift = lb.transpose() * v2;
  ASSERT_TRUE(goldShift.isApprox(s));
}
