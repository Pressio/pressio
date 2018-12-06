
#ifndef QR_FORWARD_DECLARATIONS_HPP_
#define QR_FORWARD_DECLARATIONS_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_algorithms_tags.hpp"

namespace rompp{  namespace qr{

template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type,
	 typename algorithm = ::rompp::qr::Hacked,
	 typename enable = void>
class QRSolver;


}}//end namespace rompp::qr
#endif
