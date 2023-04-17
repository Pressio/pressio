
#include <gtest/gtest.h>

#include "pressio/type_traits.hpp"

struct FakeStateType{};
struct FakeRhsType{};
struct FakeResultOfApplyMassMatType{};
struct FakeManifoldJacType{};

struct FakeTimeType{
  operator double(){ return double{}; }

  bool operator==(const FakeTimeType&) const{ return true;}
  bool operator<(const FakeTimeType&) const{ return true;}
  bool operator<=(const FakeTimeType&) const{ return true;}
  bool operator>(const FakeTimeType&) const{ return true;}
  bool operator>=(const FakeTimeType&) const{ return true;}
  bool operator!=(const FakeTimeType&) const{ return true;}
};

namespace pressio{
template<> struct Traits<FakeStateType>{
  using scalar_type = double;
  static constexpr int rank = 1;
};
template<> struct Traits<FakeRhsType>{
  using scalar_type = double;
  static constexpr int rank = 1;
};
template<> struct Traits<FakeResultOfApplyMassMatType>{
  using scalar_type = double;
  static constexpr int rank = 2;
};
template<> struct Traits<FakeManifoldJacType>{
  using scalar_type = double;
  static constexpr int rank = 2;
};
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

  FakeResultOfApplyMassMatType createResultOfMassMatrixActionOn(const FakeManifoldJacType & ) const;

  void applyMassMatrix(const state_type &,
		     const FakeManifoldJacType &,
		     const time_type &,
		     FakeResultOfApplyMassMatType &) const;
};

struct System1{
  using time_type = FakeTimeType;
  using state_type = FakeStateType;
  using right_hand_side_type = FakeRhsType;

  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type & /*unused*/,
		     const time_type & /*unused*/,
		     right_hand_side_type & /*unused*/) const{}

  FakeResultOfApplyMassMatType createResultOfMassMatrixActionOn(const FakeManifoldJacType & ) const;

  void applyMassMatrix(const state_type &,
		     const FakeManifoldJacType &,
		     const time_type &,
		     FakeResultOfApplyMassMatType &) const;
};

struct System2{
  using time_type = FakeTimeType;
  using state_type = FakeStateType;
  using right_hand_side_type = FakeRhsType;

  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type & /*unused*/,
		     const time_type & /*unused*/,
		     right_hand_side_type & /*unused*/) const{}

  FakeResultOfApplyMassMatType createResultOfMassMatrixActionOn(const FakeManifoldJacType & ) const;

  // invalidate by commenting this out
  // void applyMassMatrix(const state_type &,
  // 		     const FakeManifoldJacType &,
  // 		     const time_type &,
  // 		     FakeResultOfApplyMassMatType &) const;
};

TEST(rom, concept_semidiscrete_fom_with_jacobian_action)
{
  using namespace pressio::rom;

#ifdef PRESSIO_ENABLE_CXX20
  static_assert( SemiDiscreteFomWithMassMatrixAction<System0, FakeManifoldJacType>, "");
  static_assert( SemiDiscreteFomWithMassMatrixAction<System1, FakeManifoldJacType>, "");
  static_assert(!SemiDiscreteFomWithMassMatrixAction<System2, FakeManifoldJacType>, "");
#else
  static_assert( SemiDiscreteFomWithMassMatrixAction<System0, FakeManifoldJacType>::value, "");
  static_assert( SemiDiscreteFomWithMassMatrixAction<System1, FakeManifoldJacType>::value, "");
  static_assert(!SemiDiscreteFomWithMassMatrixAction<System2, FakeManifoldJacType>::value, "");
#endif
}
