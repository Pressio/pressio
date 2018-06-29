
#ifndef ODE_EULER_IMPLICIT_JACOBIAN_IMPL_HPP_
#define ODE_EULER_IMPLICIT_JACOBIAN_IMPL_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace impl{

// //******************************************************
// // enable for SPARSE matrix
// //******************************************************
// template <typename jacobian_type,
// 	  typename time_type,
// 	  typename
// 	  std::enable_if<core::details::traits<jacobian_type>::isMatrix==1 &&
// 			 core::details::traits<jacobian_type>::isSparse==1
// 			 >::type * = nullptr
// 	  >
// void implicit_euler_jacobian_impl(jacobian_type & jac,
// 				  time_type dt)
// {
//   using o_t =
//     typename core::details::traits<jacobian_type>::ordinal_t;
//   using sc_t =
//     typename core::details::traits<jacobian_type>::scalar_t;

//   // do first entry
//   o_t j = 0;
//   sc_t val1 = jac(0,0);
//   jac.insertValues(0, 1, &val1, &j);

//   // do others
//   o_t js[2];
//   sc_t vals[2];
//   for (size_t i=1; i < jac.rows(); ++i){
//     js[0] = i-1;
//     vals[0] = - dt * jac(i,i-1);
//     js[1] = i;
//     vals[1] = 1.0 - dt * jac(i,i);
//   }
//   // // obviously this needs to be fixed to use operator []
//   // // auto & jac = J.getNonConstRefToData();
//   // jac[0,0] = 1.0 - dt * jac[0,0];
//   // for (size_t i=1; i < jac.rows(); ++i){
//   //   jac[i,i-1] = - dt * jac[i,i-1];
//   //   jac[i,i] = 1.0 - dt * jac[i,i];
//   // }
// }

  
//******************************************************
// enable for DENSE jacobian (even though this might never
// be used because usually jacobians are sparse
//******************************************************
template <typename jacobian_type,
	  typename time_type,
	  typename
	  std::enable_if<core::details::traits<jacobian_type>::isMatrix==1 &&
			 core::details::traits<jacobian_type>::isDense==1
			 >::type * = nullptr
	  >
void implicit_euler_jacobian_impl(jacobian_type & jac,
				  time_type dt)
{
  jac(0,0) = 1.0 - dt * jac(0,0);
  for (decltype(jac.rows()) i=1; i < jac.rows(); ++i){
    jac(i,i-1) = - dt * jac(i,i-1);
    jac(i,i) = 1.0 - dt * jac(i,i);
  }
}

  
}//end namespace impl
}//end namespace ode
#endif 

