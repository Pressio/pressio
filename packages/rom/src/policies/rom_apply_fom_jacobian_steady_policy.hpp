
#ifndef ROM_APPLY_FOM_JACOBIAN_STEADY_HPP_
#define ROM_APPLY_FOM_JACOBIAN_STEADY_HPP_

namespace pressio{ namespace rom{ namespace policy{

template <>
struct ApplyFomJacobianDefault<true>{

  template <
    typename fom_t,
    typename state_t,
    typename operand_t
    >
  auto evaluate(const fom_t	  & fomObj,
		const state_t   & yFOM,
		const operand_t & B) const
    -> decltype(fomObj.applyJacobian(*yFOM.data(), *B.data()))
  {
    return fomObj.applyJacobian(*yFOM.data(), *B.data());
  }

  template <
    typename fom_t,
    typename state_t,
    typename operand_t,
    typename result_t
    >
  void evaluate(const fom_t	  & fomObj,
		const state_t	  & yFOM,
		const operand_t & B,
		result_t	  & out) const
  {
    fomObj.applyJacobian(*yFOM.data(), *B.data(), *out.data());
  }

};

}}} //end namespace pressio::rom::policy
#endif
