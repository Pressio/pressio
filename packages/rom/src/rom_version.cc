
#include "rom_version.hpp"

namespace pressio{ namespace rom{

inline std::string version(){ 
	return("rom in PRESSIO " PRESSIO_VERSION_STRING); 
}

}}//end namespace 
