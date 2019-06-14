
#ifndef ROM_EVALUATE_FOM_RHS_UNSTEADY_HPP_
#define ROM_EVALUATE_FOM_RHS_UNSTEADY_HPP_

namespace rompp{ namespace rom{ namespace policy{

template <>
struct EvaluateFomRhsDefault<false>{

  //------------------------------------------
  // enabled for native c++
  //------------------------------------------
  template <
    typename fom_t,   typename state_t,
    typename rhs_t, typename time_t
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      mpl::not_same<fom_t, pybind11::object>::value and
      !::rompp::core::meta::is_array_pybind11<state_t>::value and
      !::rompp::core::meta::is_array_pybind11<rhs_t>::value
      > * = nullptr
#endif
    >
  void evaluate(const fom_t	& fomObj,
		const state_t & yFOM,
		rhs_t		& rhs,
		time_t		t) const{
    fomObj.residual(*yFOM.data(), *rhs.data(), t);
  }

  template <
    typename fom_t, typename state_t, typename time_t
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
	mpl::not_same<fom_t, pybind11::object>::value and
	!::rompp::core::meta::is_array_pybind11<state_t>::value
	> * = nullptr
#endif
    >
  auto evaluate(const fom_t	& fomObj,
		const state_t & yFOM,
		time_t		t) const
    -> decltype(fomObj.residual(*yFOM.data(), t))
  {
    return fomObj.residual(*yFOM.data(), t);
  }


  //------------------------------------------
  // enabled when interfacing with python
  //------------------------------------------
#ifdef HAVE_PYBIND11
  template <
    typename fom_t, typename state_t,
    typename rhs_t, typename time_t,
    mpl::enable_if_t<
      mpl::is_same<fom_t, pybind11::object>::value and
      ::rompp::core::meta::is_cstyle_array_pybind11<state_t>::value and
      ::rompp::core::meta::is_cstyle_array_pybind11<rhs_t>::value and
      // because we should have all = pybind11::array_t
      mpl::is_same<state_t, rhs_t>::value
      > * = nullptr
    >
  void evaluate(const fom_t	& fomObj,
		const state_t & yFOM,
		rhs_t		& rhs,
		time_t		t) const{
    fomObj.attr("residual2")(yFOM, rhs, t);
  }

  template <
    typename fom_t, typename state_t, typename time_t,
    mpl::enable_if_t<
      mpl::is_same<fom_t, pybind11::object>::value and
      ::rompp::core::meta::is_cstyle_array_pybind11<state_t>::value
      > * = nullptr
    >
  state_t evaluate(const fom_t	& fomObj,
		   const state_t & yFOM,
		   time_t t) const {
    return fomObj.attr("residual1")(yFOM, t);
  }
#endif
};

}}} //end namespace rompp::rom::policy
#endif
