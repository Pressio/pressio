
#ifndef rom_query_fom_discrete_time_residual_HPP_
#define rom_query_fom_discrete_time_residual_HPP_

namespace pressio{ namespace rom{

// template <typename fom_state_t, typename fom_t>
// auto queryFomTimeDiscreteResidual(const fom_state_t & fomCurrentState,
// 				  const fom_t   & fomObj)
//   -> decltype(fomObj.createTimeDiscreteResidualObject(*fomCurrentState.data()))
// {
//   return fomObj.createTimeDiscreteResidualObject(*fomCurrentState.data());
// }

template <
  typename fom_state_t, typename fom_t, typename step_t, 
  typename time_t, typename norm_t, typename result_t
  >
void queryFomDiscreteTimeResidual(const fom_state_t & state_n,
				  const fom_state_t & state_nm1,
				  const fom_t   & fomObj,
				  const time_t  & time,
				  const time_t  & dt,
				  const step_t  & step,
				  result_t & R,
				  const pressio::Norm & normKind,
				  norm_t & normValue)
{
  fomObj.template discreteTimeResidual(step, time, dt, *R.data(), 
  	normKind, normValue, *state_n.data(), *state_nm1.data());
}

template <
  typename fom_state_t, typename fom_t, typename step_t, 
  typename time_t, typename norm_t, typename result_t
  >
void queryFomDiscreteTimeResidual(const fom_state_t & state_n,
				  const fom_state_t & state_nm1,
				  const fom_state_t & state_nm2,
				  const fom_t   & fomObj,
				  const time_t  & time,
				  const time_t  & dt,
				  const step_t  & step,
				  result_t & R,
				  const pressio::Norm & normKind,
				  norm_t & normValue)
{
  fomObj.template discreteTimeResidual(step, time, dt, *R.data(),
  	normKind, normValue, *state_n.data(), *state_nm1.data(), *state_nm2.data());
}

}} //end namespace pressio::rom
#endif
