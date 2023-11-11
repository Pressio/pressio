
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/rom_subspaces.hpp"

namespace impl{

template <
  class FomLinearOperatorType,
  class LeftBasisType,
  class RightBasisType
  >
class RomLinearThingy{
  using fomop_type = std::remove_cvref_t<FomLinearOperatorType>;
  using lb_type    = std::remove_cvref_t<LeftBasisType>;
  using rb_type    = std::remove_cvref_t<RightBasisType>;
  using fom_op_apply_result_type = typename FomLinearOperatorType::action_type;

  const fomop_type * fop_ = nullptr;
  const lb_type * lb_ = nullptr;
  const rb_type * rb_ = nullptr;
  fom_op_apply_result_type fomOpAction_;

  explicit RomLinearThingy(const fomop_type & fop,
			   const lb_type & lb,
			   const rb_type & rb)
    : fop_(&fop), lb_(&lb), rb_(&rb),
      fomOpAction_(fop.create){}

  void compute(phiLeft, phiRight, xRef)
  {
    const std::size_t nCols = pressio::ops::extent(*rb_);
    for (int j=0; i<nCols; ++){
      auto rb_j = pressio::ops::col(*rb_);
      fop_->apply(rb_j, fomAction);
    }

    // allocate mem for M = A*phiRight
    // loop j over col of phiRight
    //   compute Action of A on that col and store into M[:,j]
    // result = phiLeft^T M
  }

  // void export(file, ...){}
  // ... getReducedOperator(){ return result; }
  // ... getReducedOffset(){ return }
  // ... getReducedForcing(){ return ; }
};
}

template <
  class FomLinearOperatorType,
  class LeftBasisType,
  class RightBasisType
  >
auto create_linear_rom(const FomLinearOperatorType & fomOp,
		       const LeftBasisType & lb,
		       const RightBasisType & rb)
{
  using ret_t = impl::RomLinearThingy<
    FomLinearOperatorType, LeftBasisType, RightBasisType>;
  return ret_t(fomOp, ls, rb);
}

// ---------------------------------------------

class MyFomLinearOperator{
  using action_result_type = Eigen::VectorXd;

  template<class OperandT>
  action_result_type createResultOfAction(const OperandT &) const{}

  template<class OperandT>
  void apply(const OperandT & operand,
	     action_result_type & result)
  {
    // compute result = A*x
  }
};

TEST(rom, linear_rom)
{
  // using basis_t = Eigen::MatrixXd;
  // using space_t = LinearSubspace<basis_t>;
  // basis_t A = basis_t::Random(m, n);
  // space_t space(A, space_t::SpanningSet::Columns);

  // const auto & B = space.basis();
  // EXPECT_TRUE(B.data() != A.data());
  // EXPECT_TRUE( B.isApprox(A) );
  // EXPECT_TRUE( space.dimension() ==n );
  // EXPECT_TRUE(  space.isColumnSpace() );
  //EXPECT_FALSE( space.isRowSpace() );
}
