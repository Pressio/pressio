
#include "mpl_version.hpp"

namespace rompp{ namespace mpl{

inline std::string version(){ 
	return("mpl in ROMPP " ROMPP_VERSION_STRING); 
}

}}//end namespace 
