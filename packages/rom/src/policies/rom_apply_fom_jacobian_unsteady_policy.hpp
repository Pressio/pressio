
#ifndef ROM_APPLY_FOM_JACOBIAN_UNSTEADY_HPP_
#define ROM_APPLY_FOM_JACOBIAN_UNSTEADY_HPP_

namespace rompp{ namespace rom{ namespace policy{

template <>
struct ApplyFomJacobianDefault<false>{

  //------------------------------------------
  // enabled for native c++
  //------------------------------------------
  template <
    typename fom_t, typename state_t,
    typename operand_t, typename time_t
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      mpl::not_same<fom_t, pybind11::object>::value and
      !::rompp::containers::meta::is_array_pybind11<state_t>::value and
      !::rompp::containers::meta::is_array_pybind11<operand_t>::value
      > * = nullptr
#endif
    >
  auto evaluate(const fom_t	  & fomObj,
		const state_t   & yFOM,
		const operand_t & B,
		time_t		  t) const
    -> decltype(fomObj.applyJacobian(*yFOM.data(), *B.data(), t))
  {
    return fomObj.applyJacobian(*yFOM.data(), *B.data(), t);
  }

  template <
    typename fom_t, typename state_t, typename operand_t,
    typename result_t, typename time_t
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      mpl::not_same<fom_t, pybind11::object>::value and
      !::rompp::containers::meta::is_array_pybind11<state_t>::value and
      !::rompp::containers::meta::is_array_pybind11<operand_t>::value
      > * = nullptr
#endif
    >
  void evaluate(const fom_t	  & fomObj,
		const state_t	  & yFOM,
		const operand_t & B,
		result_t	  & out,
		time_t		  t) const{
    fomObj.applyJacobian(*yFOM.data(), *B.data(), *out.data(), t);
  }


#ifdef HAVE_PYBIND11
  //------------------------------------------
  // enabled when interfacing with python
  //------------------------------------------
  template <
    typename fom_t, typename state_t,
    typename operand_t, typename time_t
    , mpl::enable_if_t<
      mpl::is_same<fom_t, pybind11::object>::value and
      ::rompp::containers::meta::is_array_pybind11<state_t>::value and
      ::rompp::containers::meta::is_array_pybind11<operand_t>::value and
      // because we should have all = pybind11::array_t
      mpl::is_same<state_t, operand_t>::value
      > * = nullptr
    >
  state_t evaluate(const fom_t	  & fomObj,
		   const state_t   & yFOM,
		   const operand_t & B,
		   time_t		  t) const{
    return fomObj.attr("applyJacobian1")(yFOM, B, t);
  }

  template <
    typename fom_t, typename state_t, typename operand_t,
    typename result_t, typename time_t
    , mpl::enable_if_t<
      mpl::is_same<fom_t, pybind11::object>::value and
      ::rompp::containers::meta::is_array_pybind11<state_t>::value and
      ::rompp::containers::meta::is_array_pybind11<operand_t>::value and
      // because we should have all = pybind11::array_t
      mpl::is_same<state_t, operand_t>::value
      > * = nullptr
    >
  void evaluate(const fom_t	  & fomObj,
		const state_t	  & yFOM,
		const operand_t & B,
		result_t	  & out,
		time_t		  t) const{
    fomObj.attr("applyJacobian2")(yFOM, B, out, t);
  }
#endif
  
};

}}} //end namespace rompp::rom::policy
#endif
