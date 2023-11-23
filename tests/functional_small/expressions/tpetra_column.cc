
#include <gtest/gtest.h>
#include "pressio/expressions.hpp"
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorGlobSize15Fixture, multi_vector_column_expr)
{
  {
    auto expr = pressio::column(*myMv_, 0);
    auto v = expr.native();
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());

    std::vector<ST> gold(localSize_);
    if (rank_ == 0){
      gold = {1,5,9,13,17};
    }
    else if (rank_==1){
      gold = {21,25,29,33,37};
    }
    else{
      gold = {41,45,49,53,57};
    }

    for (std::size_t i=0; i<expr.extentLocal(0); ++i){
      ASSERT_TRUE(v_h(i,0) == gold[i]);
    }
  }

  {
    auto expr = pressio::column(*myMv_, 3);
    auto v = expr.native();
    auto v_h = v.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());

    std::vector<ST> gold(localSize_);
    if (rank_ == 0){
      gold = {4,8,12,16,20};
    }
    else if (rank_==1){
      gold = {24,28,32,36,40};
    }
    else{
      gold = {44,48,52,56,60};
    }

    for (std::size_t i=0; i<expr.extentLocal(0); ++i){
      ASSERT_TRUE(v_h(i,0) == gold[i]);
    }
  }

}
