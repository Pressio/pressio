
#include <gtest/gtest.h>
#include "ode_meta.hpp"


TEST(ode_meta, checkStateType)
{
  using st_t1 = std::vector<double>;
  using st_t2 = double;
  using st_t3 = std::map<int,int>;

  static_assert( ode::meta::isLegitimateStateType<st_t1>::value,"");
  static_assert( !ode::meta::isLegitimateStateType<st_t2>::value,"");
  static_assert( ode::meta::isLegitimateStateType<st_t3>::value,"");
}
//--------------------------------------------------------------------

TEST(ode_meta, checkResidualType)
{
  using st_t1 = std::vector<double>;
  using st_t2 = double;
  using st_t3 = std::map<int,int>;

  static_assert( ode::meta::isLegitimateResidualType<st_t1>::value,"");
  static_assert( !ode::meta::isLegitimateResidualType<st_t2>::value,"");
  static_assert( ode::meta::isLegitimateResidualType<st_t3>::value,"");  
}
//--------------------------------------------------------------------

TEST(ode_meta, checkTimeType)
{
  using st_t1 = std::vector<double>;
  using st_t2 = double;
  using st_t3 = std::map<int,int>;
  using st_t4 = int;

  static_assert( !ode::meta::isLegitimateTimeType<st_t1>::value,"");
  static_assert( ode::meta::isLegitimateTimeType<st_t2>::value,"");
  static_assert( !ode::meta::isLegitimateTimeType<st_t3>::value,"");  
  static_assert( !ode::meta::isLegitimateTimeType<st_t4>::value,"");
}
//--------------------------------------------------------------------
