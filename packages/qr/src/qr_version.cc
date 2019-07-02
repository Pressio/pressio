
#include "qr_version.hpp"

namespace pressio{ namespace qr{

inline std::string version(){ 
	return("qr in PRESSIO " PRESSIO_VERSION_STRING); 
}

}}//end namespace 
