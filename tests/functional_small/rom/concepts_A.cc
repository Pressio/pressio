
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

struct FakeStateType{};
struct FakeRhsType{};
struct FakeTimeType{
  operator double(){ return double{}; }
};

namespace pressio{
template<> struct Traits<FakeStateType>{ using scalar_type = double; };
template<> struct Traits<FakeRhsType>{ using scalar_type = double; };
}

#include "pressio/rom_concepts.hpp"

//
// rhs only
//
struct System0{
  using time_type = double;
  using state_type = FakeStateType;
  using right_hand_side_type = FakeRhsType;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type & /*unused*/,
		     const double & /*unused*/,
		     right_hand_side_type & /*unused*/) const{}
};

struct System1{
  using time_type = FakeTimeType;
  using state_type = FakeStateType;
  using right_hand_side_type = FakeRhsType;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type & /*unused*/,
		     time_type /*unused*/,
		     right_hand_side_type & /*unused*/) const{}
};

struct System2{
  using time_type = FakeTimeType;
  using state_type = FakeStateType;
  using right_hand_side_type = FakeRhsType;

  // invalidate by omitting createState
  //state_type createState() const{ return state_type(); }

  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type & /*unused*/,
		     time_type /*unused*/,
		     right_hand_side_type & /*unused*/) const{}
};

TEST(rom, concepts_A)
{
  using namespace pressio::rom;

  static_assert( SemiDiscreteFom<System0>::value, "");
  static_assert( SemiDiscreteFom<System1>::value, "");
  static_assert(!SemiDiscreteFom<System2>::value, "");
}
