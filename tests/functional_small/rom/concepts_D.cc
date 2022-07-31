
#include <gtest/gtest.h>
#include "pressio/rom_concepts.hpp"

struct FakeType1{};
struct FakeType2{};

#define NT1() using operand_type = FakeType1;
#define NT2() using result_type  = FakeType2;

#define M1() result_type createApplyMaskResult(const operand_type &) const;
#define M2() void operator()(const operand_type &, result_type &) const;

struct S1{
  NT1() NT2() M1() M2()
};

struct S2{
  NT1() NT2() M2()
};

struct S3{
  NT1() NT2() M1()
};

TEST(rom, concepts_D)
{
  using namespace pressio::rom;
  static_assert( TimeInvariantMasker<S1>::value, "");
  static_assert( !TimeInvariantMasker<S2>::value, "");
  static_assert( !TrialSubspace<S3>::value, "");
}
