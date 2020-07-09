
#ifndef ROM_QUERY_FOM_APPLY_DISCRETE_TIME_JACOBIAN_HPP_
#define ROM_QUERY_FOM_APPLY_DISCRETE_TIME_JACOBIAN_HPP_

namespace pressio{ namespace rom{

// template <class fom_state_t, class fom_t, class operand_t>
// auto queryFomApplyTimeDiscreteJacobian(const fom_state_t & fomCurrentState,
// 				       const fom_t & fomObj,
// 				       const operand_t & B)
//   -> decltype(fomObj.createApplyTimeDiscreteJacobianObject(*fomCurrentState.data(), *B.data()))
// {
//   return fomObj.createApplyTimeDiscreteJacobianObject(*fomCurrentState.data(),
// 						      *B.data());
// }

template <
  class fom_state_t, class fom_t, class step_t, class time_t, class operand_t, class result_t
  >
void queryFomApplyDiscreteTimeJacobian(const fom_state_t & state_n,
				       const fom_state_t & state_nm1,
				       const fom_t & fomObj,
				       const time_t & time,
				       const time_t & dt,
				       const step_t & step,
				       const operand_t & operand,
				       result_t & result)
{
  fomObj.template applyDiscreteTimeJacobian(step, time, dt,
					    *operand.data(),
					    *result.data(),
					    *state_n.data(),
					    *state_nm1.data());
}

template <
  class fom_state_t, class fom_t, class step_t, class time_t, class operand_t, class result_t
  >
void queryFomApplyDiscreteTimeJacobian(const fom_state_t & state_n,
				       const fom_state_t & state_nm1,
				       const fom_state_t & state_nm2,
				       const fom_t & fomObj,
				       const time_t & time,
				       const time_t & dt,
				       const step_t & step,
				       const operand_t & operand,
				       result_t & result)
{
  fomObj.template applyDiscreteTimeJacobian(step, time, dt,
					    *operand.data(),
					    *result.data(),
					    *state_n.data(),
					    *state_nm1.data(),
					    *state_nm2.data());
}

}} //end namespace pressio::rom
#endif
