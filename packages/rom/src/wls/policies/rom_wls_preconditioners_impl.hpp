#ifndef ROM_WLS_PRECONDITIONERS_IMPL_HPP_
#define ROM_WLS_PRECONDITIONERS_IMPL_HPP_
/*
Preconditioner objects for WLS. These act on the time local residuals and Jacobians
*/
namespace pressio{ namespace rom{ namespace wls{ namespace preconditioners{ namespace impl{

struct NoPreconditioner{
template <typename app_t, 
          typename fom_state_t, 
          typename operand_t, 
          typename scalar_t>
void operator()(const app_t &appObj,
                const fom_state_t & yFom, 
                operand_t & operand,
                const scalar_t & t)const{}
};


struct AppPreconditioner{

template <typename app_t, 
          typename fom_state_t, 
          typename operand_t, 
          typename scalar_t>
void operator()(const app_t & appObj, 
                const fom_state_t & yFom, 
                operand_t & operand,
                const scalar_t &t) const
{
    appObj.applyPreconditioner(*yFom.data(), *operand.data(), t);
}
};
}}}}}//end namespace pressio
#endif
