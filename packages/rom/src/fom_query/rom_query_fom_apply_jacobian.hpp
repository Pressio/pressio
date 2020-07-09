
#ifndef ROM_QUERY_FOM_APPLY_JACOBIAN_HPP_
#define ROM_QUERY_FOM_APPLY_JACOBIAN_HPP_

namespace pressio{ namespace rom{

template <
  typename fom_t, typename state_t, typename operand_t, typename result_t, typename time_type
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  , mpl::enable_if_t<
      mpl::not_same<fom_t, pybind11::object>::value and
      !::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value and
      !::pressio::containers::predicates::is_matrix_wrapper_pybind<operand_t>::value,
      int > = 0
#endif
  >
void queryFomApplyJacobian(const fom_t & fomObj,
           const state_t & fomState,
           const operand_t & operand,
           result_t & result,
           const time_type & time)
{
  fomObj.applyJacobian(*fomState.data(), *operand.data(), time, *result.data());
}

template <
  typename fom_t, typename state_t, typename operand_t, typename result_t
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  , mpl::enable_if_t<
      mpl::not_same<fom_t, pybind11::object>::value and
      !::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value and
      !::pressio::containers::predicates::is_matrix_wrapper_pybind<operand_t>::value,
      int > = 0
#endif
  >
void queryFomApplyJacobian(const fom_t & fomObj,
           const state_t & fomState,
           const operand_t & operand,
           result_t & result)
{
  fomObj.applyJacobian(*fomState.data(), *operand.data(), *result.data());
}



// template <
//   typename fom_t, typename state_t, typename operand_t, typename time_t
// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//   , mpl::enable_if_t<
//       mpl::not_same<fom_t, pybind11::object>::value and
//       !::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value and
//       !::pressio::containers::predicates::is_matrix_wrapper_pybind<operand_t>::value,
//       int > = 0
// #endif
//   >
// auto queryFomApplyJacobianUnsteady(const fom_t & fomObj,
// 				   const state_t & yFOM,
// 				   const operand_t & B,
// 				   const time_t & t)
//   -> decltype(fomObj.applyJacobian(*yFOM.data(), *B.data(), t))
// {
//   return fomObj.applyJacobian(*yFOM.data(), *B.data(), t);
// }


// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template <typename state_t, typename operand_t, typename result_t, typename time_t>
// mpl::enable_if_t<
//   ::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value and
//   ::pressio::containers::predicates::is_matrix_wrapper_pybind<operand_t>::value
//   >
// queryFomApplyJacobianUnsteady(const pybind11::object  & fomObj,
// 			      const state_t	  & yFOM,
// 			      const operand_t & B,
// 			      result_t	  & out,
// 			      const time_t & t)
// {
//   *out.data() = fomObj.attr("applyJacobian")(*yFOM.data(), *B.data(), t);
// }

// template <typename state_t, typename operand_t, typename time_t>
// mpl::enable_if_t<
//   ::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value and
//   ::pressio::containers::predicates::is_matrix_wrapper_pybind<operand_t>::value,
//   typename ::pressio::containers::details::traits<state_t>::wrapped_t
//   >
// queryFomApplyJacobianUnsteady(const pybind11::object & fomObj,
// 			      const state_t   & yFOM,
// 			      const operand_t & B,
// 			      const time_t & t)
// {
//   return fomObj.attr("applyJacobian")(*yFOM.data(), *B.data(), t);
// }
// #endif

}} //end namespace pressio::rom
#endif
