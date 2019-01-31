
#ifndef ROM_APPLY_FOM_JACOBIAN_HPP_
#define ROM_APPLY_FOM_JACOBIAN_HPP_

namespace rompp{ namespace rom{ namespace policy{

struct ApplyFomJacobianDefault{

  template <
    typename fom_t,
    typename state_w_t,
    typename operand_w_t,
    typename time_t
    >
  auto evaluate(const fom_t	  & fomObj,
		const state_w_t   & yFOM,
		const operand_w_t & B,
		time_t		  t) const
    -> decltype(fomObj.applyJacobian(*yFOM.data(), *B.data(), t))
  {
    return fomObj.applyJacobian(*yFOM.data(), *B.data(), t);
  }

  template <
    typename fom_t,
    typename state_w_t,
    typename operand_w_t,
    typename result_t,
    typename time_t
    >
  void evaluate(const fom_t	  & fomObj,
		const state_w_t	  & yFOM,
		const operand_w_t & B,
		result_t	  & out,
		time_t		  t) const
  {
    fomObj.applyJacobian(*yFOM.data(), *B.data(), *out.data(), t);
  }

};

}}} //end namespace rompp::rom::policy
#endif
