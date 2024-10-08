
#ifndef SOLVERS_NONLINEAR_IMPL_INTERNAL_TAGS_HPP_
#define SOLVERS_NONLINEAR_IMPL_INTERNAL_TAGS_HPP_

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

// this tag is inside the impl namespace because we do NOT want
// to expose it outside, this is an impl detail
struct SystemTag{};

// this represent ths Q^T *r for QR gauss newton
struct QTransposeResidualTag{};


struct NewtonTag{};
struct GaussNewtonNormalEqTag{};
struct WeightedGaussNewtonNormalEqTag{};
struct LevenbergMarquardtNormalEqTag{};
struct GaussNewtonQrTag{};

}}}
#endif
