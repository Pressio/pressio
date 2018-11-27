
#ifndef QR_FORWARD_DECLARATIONS_HPP_
#define QR_FORWARD_DECLARATIONS_HPP_

#include "qr_ConfigDefs.hpp"

namespace rompp{  namespace qr{ namespace hack{

template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type,
	 typename enable = void>
class QRSolver;


}}}//end namespace rompp::qr::hack
#endif
