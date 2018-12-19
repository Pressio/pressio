
#ifndef QR_FORWARD_DECLARATIONS_HPP_
#define QR_FORWARD_DECLARATIONS_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_algorithms_tags.hpp"
#include "../../core/src/core_forward_declarations.hpp"

namespace rompp{  namespace qr{

template<typename matrix_type,
	 typename R_type,
	 typename algorithm,
	 template <typename...> class Q_type =
	 ::rompp::core::MultiVector,
	 typename enable = void>
class QRSolver;


}}//end namespace rompp::qr
#endif
