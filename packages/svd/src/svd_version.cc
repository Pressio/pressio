
#include "svd_version.hpp"

namespace rompp{ namespace svd{

inline std::string version(){ 
	return("svd in ROMPP " ROMPP_VERSION_STRING); 
}

}}//end namespace 
