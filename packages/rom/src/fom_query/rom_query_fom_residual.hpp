
#ifndef ROM_QUERY_FOM_RESIDUAL_HPP_
#define ROM_QUERY_FOM_RESIDUAL_HPP_

namespace pressio{ namespace rom{

template <typename fom_t, typename state_t, typename result_t>
void queryFomResidual(const fom_t & fomObj,
			    const state_t & yFOM,
			    result_t & R)
{
  fomObj.residual(*yFOM.data(), *R.data());
}

// template <typename fom_t, typename state_t>
// auto queryFomVelocitySteady(const fom_t & fomObj,
// 			    const state_t & yFOM)
//   -> decltype(fomObj.velocity(*yFOM.data()))
// {
//   return fomObj.velocity(*yFOM.data());
// }

}} //end namespace pressio::rom
#endif
