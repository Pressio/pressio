
#ifndef QR_ALGORITHMS_TAGS_HPP_
#define QR_ALGORITHMS_TAGS_HPP_

namespace rompp{ namespace qr{

struct GramSchmidt{};
struct Householder{};

#if defined HAVE_TRILINOS and defined HAVE_ANASAZI_TSQR
//this defaults to Anasazi
struct TSQR{};
#endif

#if defined HAVE_TRILINOS and defined HAVE_BELOS_TSQR
struct TSQRBelos{};
#endif


}} // end namespace rompp::qr
#endif
