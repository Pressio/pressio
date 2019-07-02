
#include "svd_version.hpp"

namespace pressio{ namespace svd{

inline std::string version(){ 
	return("svd in PRESSIO " PRESSIO_VERSION_STRING); 
}

}}//end namespace 
