
#include "rom_version.hpp"

namespace rompp{ namespace rom{

inline std::string version(){ 
	return("rom in ROMPP " ROMPP_VERSION_STRING); 
}

}}//end namespace 
