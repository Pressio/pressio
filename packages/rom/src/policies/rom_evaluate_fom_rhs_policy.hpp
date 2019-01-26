
#ifndef ROM_EVALUATE_FOM_RHS_HPP_
#define ROM_EVALUATE_FOM_RHS_HPP_

namespace rompp{ namespace rom{ namespace policy{

struct EvaluateFomRhsDefault{

  template <
    typename app_t,
    typename fom_st_t,
    typename result_t,
    typename time_t
    >
  void evaluate(const app_t    & app,
		const fom_st_t & yFOM,
		result_t       & rhs,
		time_t	       t) const
  {
    app.residual(*yFOM.data(), *rhs.data(), t);
  }


  template <
    typename app_t,
    typename fom_st_t,
    typename time_t
    >
  auto evaluate(const app_t    & app,
		const fom_st_t & yFOM,
		time_t	       t) const
    -> decltype(app.residual(*yFOM.data(), t))
  {
    return app.residual(*yFOM.data(), t);
  }

};

}}} //end namespace rompp::rom::policy
#endif
