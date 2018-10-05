
#include <gtest/gtest.h>
#include "ode_meta.hpp"
#include "vector/core_vector_distributed_epetra.hpp"
#include "vector/core_vector_serial_eigen.hpp"
#include "vector/core_vector_serial_stdlib.hpp"

TEST(ode_meta, checkStateType)
{
  using st_t1 = std::vector<double>;
  using st_t2 = double;
  using st_t3 = std::map<int,int>;
  static_assert( core::meta::is_vector_stdlib<st_t1>::value, "");
  static_assert( ode::meta::is_legitimate_explicit_state_type<st_t1>::value,"");
  using myt = core::Vector<st_t1>;
  static_assert( ode::meta::is_legitimate_explicit_state_type<myt>::value,"");  
  static_assert( !ode::meta::is_legitimate_explicit_state_type<st_t2>::value,"");
  static_assert( !ode::meta::is_legitimate_explicit_state_type<st_t3>::value,"");
}
//--------------------------------------------------------------------

TEST(ode_meta, checkResidualType)
{
  using st_t1 = std::vector<double>;
  using st_t2 = double;
  using st_t3 = std::map<int,int>;
  static_assert( ode::meta::is_legitimate_explicit_residual_type<st_t1>::value,"");
  static_assert( !ode::meta::is_legitimate_explicit_residual_type<st_t2>::value,"");
  static_assert( !ode::meta::is_legitimate_explicit_residual_type<st_t3>::value,"");  
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
