
#include "qr_version.hpp"

namespace rompp{ namespace qr{

inline std::string version(){ 
	return("qr in ROMPP " ROMPP_VERSION_STRING); 
}

}}//end namespace 
