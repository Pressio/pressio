
#ifndef QR_ALGORITHMS_TAGS_HPP_
#define QR_ALGORITHMS_TAGS_HPP_

namespace rompp{ namespace qr{

struct GramSchmidt{};
struct Householder{};

// for trilinos, this defaults to Anasazi
struct TSQR{};

// this is only applicable to trilinos
struct TSQRBelos{};

}} // end namespace rompp::qr
#endif
