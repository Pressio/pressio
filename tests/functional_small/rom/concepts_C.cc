
#include <gtest/gtest.h>
#include "pressio/rom_concepts.hpp"

#define NT1() using reduced_state_type = Eigen::VectorXd;
#define NT2() using basis_type         = Eigen::MatrixXd;
#define NT3() using full_state_type    = Eigen::VectorXd;

#define M1() reduced_state_type createReducedState() const;
#define M2() full_state_type createFullState() const;
#define M3() void mapFromReducedState(const reduced_state_type &, full_state_type &) const;
#define M4() full_state_type createFullStateFromReducedState(const reduced_state_type &) const;
#define M5() const basis_type & viewBasis() const;
#define M6() const full_state_type & viewAffineOffset() const;

struct S1{
  NT1(); NT2(); NT3(); M1(); M2(); M3(); M4(); M5();
};

struct S2{
  NT1(); NT2(); NT3(); M1(); M2(); M3(); M4(); M5(); M6();
};

struct S3{
  NT1(); NT2(); NT3(); M1(); M2(); M3();
};

struct S4{
  NT1(); NT2(); NT3(); M1(); M4(); M5(); M6();
};

TEST(rom, concepts_C)
{
  using namespace pressio::rom;
  static_assert( TrialSubspace<S1>::value, "");
  static_assert( !AffineTrialSubspace<S1>::value, "");
  static_assert( TrialSubspace<S2>::value, "");
  static_assert( AffineTrialSubspace<S2>::value, "");

  static_assert( !TrialSubspace<S3>::value &&
		 !AffineTrialSubspace<S3>::value, "");
  static_assert( !TrialSubspace<S4>::value &&
		 !AffineTrialSubspace<S4>::value, "");
}
