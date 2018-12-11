
#ifndef ODE_ENUM_HPP_
#define ODE_ENUM_HPP_

namespace rompp{ namespace ode{

enum class ExplicitEnum{Undefined,
			Euler,
			RungeKutta4
};

enum class ImplicitEnum{Undefined,
			Euler,
			BDF2
};

}}//end namespace rompp::ode
#endif
