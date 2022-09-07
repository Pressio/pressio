
#include <gtest/gtest.h>
#include "pressio/rom_concepts.hpp"

struct FakeType1{};
struct FakeType2{};
struct FakeType3{};
struct FakeType4{};

#define NT1() using reduced_state_type = FakeType1;
#define NT2() using basis_type         = FakeType2;
#define NT3() using full_state_type    = FakeType3;
#define NT4() using offset_type        = FakeType4;

#define M1() reduced_state_type createReducedState() const;
#define M2() full_state_type createFullState() const;
#define M3() void mapFromReducedState(const reduced_state_type &, full_state_type &) const;
#define M4() full_state_type createFullStateFromReducedState(const reduced_state_type &) const;
#define M5() const basis_type & viewBasis() const;
#define M6() const offset_type & viewAffineOffset() const;

struct S1{
  NT1() NT2() NT3() M1() M2() M3() M4() M5()
};

struct S2{
  NT1() NT2() NT3() NT4() M1() M2() M3() M4() M5() M6()
};

struct S3{
  NT1() NT2() NT3() NT4() M1() M2() M3()
};

struct S4{
  NT1() NT2() NT3() NT4() M1() M4() M5() M6()
};

TEST(rom, concepts_C)
{
  using namespace pressio::rom;
  static_assert(  TrialColumnSubspaceConcept<S1>::value, "");
  static_assert( !AffineTrialColumnSubspaceConcept<S1>::value, "");

  static_assert( TrialColumnSubspaceConcept<S2>::value, "");
  static_assert( AffineTrialColumnSubspaceConcept<S2>::value, "");

  static_assert( !TrialColumnSubspaceConcept<S3>::value &&
		 !AffineTrialColumnSubspaceConcept<S3>::value, "");
  static_assert( !TrialColumnSubspaceConcept<S4>::value &&
		 !AffineTrialColumnSubspaceConcept<S4>::value, "");
}
