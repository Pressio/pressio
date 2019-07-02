
#include "mpl_version.hpp"

namespace pressio{ namespace mpl{

inline std::string version(){ 
	return("mpl in PRESSIO " PRESSIO_VERSION_STRING); 
}

}}//end namespace 
