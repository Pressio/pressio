
#ifndef ROM_EVALUATE_FOM_RHS_HPP_
#define ROM_EVALUATE_FOM_RHS_HPP_

namespace rompp{ namespace rom{ namespace policy{

struct EvaluateFomRhsDefault{

  template <
    typename fom_t,
    typename state_w_t,
    typename rhs_w_t,
    typename time_t
    >
  void evaluate(const fom_t	& fomObj,
		const state_w_t & yFOM,
		rhs_w_t		& rhs,
		time_t		t) const
  {
    fomObj.residual(*yFOM.data(), *rhs.data(), t);
  }


  template <
    typename fom_t,
    typename state_w_t,
    typename time_t
    >
  auto evaluate(const fom_t	& fomObj,
		const state_w_t & yFOM,
		time_t		t) const
    -> decltype(fomObj.residual(*yFOM.data(), t))
  {
    return fomObj.residual(*yFOM.data(), t);
  }

};

}}} //end namespace rompp::rom::policy
#endif
