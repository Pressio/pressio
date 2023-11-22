
#include <gtest/gtest.h>
#include "pressio/expressions.hpp"
#include "tpetra_block_only_fixtures.hpp"

using ft = tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture;

TEST_F(ft, multi_vector_column_expr)
{
  for (int k=0; k<3; ++k)
  {
    auto e = pressio::column(*myMv_, k);
    ASSERT_TRUE(e.native().getNumVectors()==1);

    auto a = e.native().getLocalViewHost(Tpetra::Access::ReadOnlyStruct());

    auto tpmv = myMv_->getMultiVectorView();
    auto e2 = pressio::column(tpmv, k);
    auto b = e2.native().getLocalViewHost(Tpetra::Access::ReadOnlyStruct());

    EXPECT_TRUE(a.extent(0) == b.extent(0));
    EXPECT_TRUE(a.extent(1) == b.extent(1));
    EXPECT_TRUE(a.extent(1) == 1);

    for (size_t i=0; i<a.extent(0); ++i){
      for (size_t j=0; j<a.extent(1); ++j){
	ASSERT_TRUE(a(i,j) == b(i,j));
      }
    }
  }
}
