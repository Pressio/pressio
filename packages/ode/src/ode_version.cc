
#include "ode_version.hpp"

namespace pressio{ namespace ode{

inline std::string version(){ 
	return("ode in PRESSIO " PRESSIO_VERSION_STRING); 
}

}}//end namespace 
