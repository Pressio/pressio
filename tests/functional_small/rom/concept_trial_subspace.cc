
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"

namespace{
struct FakeType2{};
struct FakeType3{};
}

namespace pressio{
template<> struct Traits<FakeType2>{
  static constexpr int rank = 2;
  using scalar_type = double;
};
template<> struct Traits<FakeType3>{
  static constexpr int rank = 1;
  using scalar_type = double;
};
}

namespace{
#define NT1() using reduced_state_type = Eigen::VectorXd;
#define NT2() using basis_matrix_type  = FakeType2;
#define NT3() using full_state_type    = FakeType3;

#define M1() reduced_state_type createReducedState() const;
#define M2() full_state_type createFullState() const;
#define M3() void mapFromReducedState(const reduced_state_type &, full_state_type &) const;
#define M4() full_state_type createFullStateFromReducedState(const reduced_state_type &) const;
#define M5() const FakeType2 & basis() const;
#define M6() const full_state_type & translationVector() const;
#define M7() int  dimension() const;
#define M8() bool isColumnSpace() const;
#define M9() bool isRowSpace() const;
#define M10() const FakeType2 & basisOfTranslatedSpace() const;

struct S1{
  S1& operator=(const S1 &) = delete;
  NT1() NT2() NT3() M1() M2() M3() M4() M5() M6() M7() M8() M9() M10()
};

struct S2{
  S2& operator=(const S2 &) = delete;
  NT1() NT2() NT3() M1() /*M2()*/ M3() M4() M5() M6() M7() M8() M9() M10()
};

struct S3{
  S3& operator=(const S3 &) = delete;
  NT1() /*NT2()*/ NT3() M1() M2() M3() M4() M5() M6() M7() M8() M9() M10()
};
}

TEST(rom_concepts, possibly_affine_trial_subspace1)
{
  using namespace pressio::rom;

#ifdef PRESSIO_ENABLE_CXX20
  static_assert(PossiblyAffineTrialColumnSubspace<S1>, "");
  static_assert(!PossiblyAffineTrialColumnSubspace<S2>, "");
  static_assert(!PossiblyAffineTrialColumnSubspace<S3>, "");
#else
  static_assert(PossiblyAffineTrialColumnSubspace<S1>::value, "");
  static_assert(!PossiblyAffineTrialColumnSubspace<S2>::value, "");
  static_assert(!PossiblyAffineTrialColumnSubspace<S3>::value, "");
#endif
}

TEST(rom_concepts, possibly_affine_trial_subspace2)
{
  using namespace pressio::rom;

  using phi_t = Eigen::Matrix<double, -1,-1>;
  using reduced_state_type = Eigen::VectorXd;
  phi_t phi(10, 3);
  Eigen::VectorXd shift(10);
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi, shift, false);

#ifdef PRESSIO_ENABLE_CXX20
  static_assert(PossiblyAffineTrialColumnSubspace<decltype(space)>, "");
#else
  static_assert(PossiblyAffineTrialColumnSubspace<decltype(space)>::value, "");
#endif
}
