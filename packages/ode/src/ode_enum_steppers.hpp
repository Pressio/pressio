
#ifndef ODE_ENUM_STEPPERS_HPP_
#define ODE_ENUM_STEPPERS_HPP_

namespace rompp{
namespace ode{

enum class ExplicitSteppersEnum{Undefined,
				Euler,
				RungeKutta4
};


enum class ImplicitSteppersEnum{Undefined,
				Euler,
				BDF2
};
  
} // end namespace ode
}//end namespace rompp
#endif
