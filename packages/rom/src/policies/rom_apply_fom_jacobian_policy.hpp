
#ifndef ROM_APPLY_FOM_JACOBIAN_HPP_
#define ROM_APPLY_FOM_JACOBIAN_HPP_

namespace rompp{ namespace rom{ namespace policy{

struct ApplyFomJacobianDefault{

  template <
    typename app_t,
    typename fom_state_t,
    typename operand_t,
    typename time_t
    >
  auto evaluate(const app_t	  & app,
		const fom_state_t & yFOM,
		const operand_t   & B,
		time_t		  t) const
    -> decltype(app.applyJacobian(*yFOM.data(), *B.data(), t))
  {
    return app.applyJacobian(*yFOM.data(), *B.data(), t);
  }


  template <
    typename app_t,
    typename fom_state_t,
    typename operand_t,
    typename result_t,
    typename time_t
    >
  void evaluate(const app_t	  & app,
		const fom_state_t & yFOM,
		const operand_t	  & B,
		result_t	  & out,
		time_t		  t) const
  {
    app.applyJacobian(*yFOM.data(), *B.data(), *out.data(), t);
  }

};

}}} //end namespace rompp::rom::policy
#endif
