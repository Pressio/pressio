
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

namespace{
struct FakeType1{};
}

namespace pressio{
template<>
struct Traits<FakeType1>
{
  static constexpr int rank = 2;
  using scalar_type = float;
};
}

namespace{
#define NT1() using basis_matrix_type  = FakeType1;

#define M1() const basis_matrix_type & basis() const;
#define M2() int dimension() const;
#define M3() bool isColumnSpace() const;
#define M4() bool isRowSpace() const;

class S1{
public:
  S1& operator=(const S1 &) = delete;
  NT1() M1() M2() M3() M4()
};

class S2{
public:
  // force it to be assignable whcih breaks concept
  //S1& operator=(const S1 &) = delete;
  NT1() M1() M2() M3() M4()
};

class S3{
public:
  S1& operator=(const S1 &) = delete;
  NT1() /*M1()*/ M2() M3() M4()
};

class S4{
public:
  S1& operator=(const S1 &) = delete;
  NT1() M1() /*M2()*/ M3() M4()
};

class S5{
public:
  S1& operator=(const S1 &) = delete;
  NT1() M1() M2() /*M3()*/ M4()
};

class S6{
public:
  S1& operator=(const S1 &) = delete;
  NT1() M1() M2() M3() /*M4()*/
};
}

#include "pressio/rom_concepts.hpp"
#include "pressio/rom_subspaces.hpp"

TEST(rom_concepts, vector_subspace1)
{
  using namespace pressio::rom;

  static_assert(VectorSubspace<S1>::value, "");
  static_assert(!VectorSubspace<S2>::value, "");
  static_assert(!VectorSubspace<S3>::value, "");
  static_assert(!VectorSubspace<S4>::value, "");
  static_assert(!VectorSubspace<S5>::value, "");
  static_assert(!VectorSubspace<S6>::value, "");
}

TEST(rom_concepts, vector_subspace2)
{
  using namespace pressio::rom;

  using T = LinearSubspace<FakeType1>;
  static_assert(VectorSubspace<T>::value, "");
}
