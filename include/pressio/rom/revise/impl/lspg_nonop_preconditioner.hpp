
#ifndef ROM_IMPL_LSPG_NOOP_PRECONDITIONER_HPP_
#define ROM_IMPL_LSPG_NOOP_PRECONDITIONER_HPP_

namespace pressio{ namespace rom{ namespace impl{

struct NoOpPreconditionerSteadyLspg{
  template<class StateType, class ResidualType, class JacobianType>
  void operator()(const StateType & /*state*/,
		  ResidualType & /*residual*/,
		  JacobianType & /*jacobian*/) const
  {
    // no op
  }
};

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_LSPG_HELPERS_HPP_
