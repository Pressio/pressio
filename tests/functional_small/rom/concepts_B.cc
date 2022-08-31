
#include <gtest/gtest.h>

#include "pressio/type_traits.hpp"

struct FakeStateType{};
struct FakeRhsType{};
struct FakeResultOfApplyJacType{};
struct FakeManifoldJacType{};
struct FakeTimeType{
  operator double(){ return double{}; }
};

namespace pressio{
template<> struct Traits<FakeStateType>{ using scalar_type = double; };
template<> struct Traits<FakeRhsType>{ using scalar_type = double; };
template<> struct Traits<FakeResultOfApplyJacType>{ using scalar_type = double; };
template<> struct Traits<FakeManifoldJacType>{ using scalar_type = double; };
}

#include "pressio/rom_concepts.hpp"

//
// rhs only
//
struct System0{
  using time_type = double;
  using state_type = FakeStateType;
  using right_hand_side_type = FakeRhsType;

  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type & /*unused*/,
		     const double & /*unused*/,
		     right_hand_side_type & /*unused*/) const{}

  FakeResultOfApplyJacType createApplyJacobianResult(const FakeManifoldJacType & ) const;

  void applyJacobian(const state_type &,
		     const FakeManifoldJacType &,
		     const time_type &,
		     FakeResultOfApplyJacType &) const;
};

struct System1{
  using time_type = FakeTimeType;
  using state_type = FakeStateType;
  using right_hand_side_type = FakeRhsType;

  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type & /*unused*/,
		     const time_type & /*unused*/,
		     right_hand_side_type & /*unused*/) const{}

  FakeResultOfApplyJacType createApplyJacobianResult(const FakeManifoldJacType & ) const;

  void applyJacobian(const state_type &,
		     const FakeManifoldJacType &,
		     const time_type &,
		     FakeResultOfApplyJacType &) const;
};

struct System2{
  using time_type = FakeTimeType;
  using state_type = FakeStateType;
  using right_hand_side_type = FakeRhsType;

  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type & /*unused*/,
		     const time_type & /*unused*/,
		     right_hand_side_type & /*unused*/) const{}

  FakeResultOfApplyJacType createApplyJacobianResult(const FakeManifoldJacType & ) const;

  // invalidate by commenting this out
  // void applyJacobian(const state_type &,
  // 		     const FakeManifoldJacType &,
  // 		     const time_type &,
  // 		     FakeResultOfApplyJacType &) const;
};

TEST(rom, concepts_B)
{
  using namespace pressio::rom;

  static_assert( SemiDiscreteFomWithJacobianAction<System0, FakeManifoldJacType>::value, "");
  static_assert( SemiDiscreteFomWithJacobianAction<System1, FakeManifoldJacType>::value, "");
  static_assert(!SemiDiscreteFomWithJacobianAction<System2, FakeManifoldJacType>::value, "");
}
