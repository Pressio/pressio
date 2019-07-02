
#include "containers_version.hpp"

namespace pressio{ namespace containers{

inline std::string version(){ 
	return("containers in PRESSIO " PRESSIO_VERSION_STRING); 
}

}}//end namespace pressio::containers
