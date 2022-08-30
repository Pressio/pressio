
#include <gtest/gtest.h>
#include "pressio_rom.hpp"

using fom_velo_t = ::pressio::containers::Vector<Eigen::VectorXd>;
using dmat_t = ::pressio::containers::DenseMatrix<Eigen::MatrixXd>;
using rom_state_t = ::pressio::containers::Vector<Eigen::VectorXd>;
using rom_jac_t = ::pressio::containers::DenseMatrix<Eigen::MatrixXd>;

struct Projector
{
  void apply(const fom_velo_t &, rom_state_t &) const{};
  void apply(const dmat_t &, rom_jac_t &) const{};
};

TEST(rom_galerkin_meta, validProjectorExplicitStepping)
{
  static_assert
    (pressio::rom::galerkin::constraints::projector_explicit_stepping
     <Projector, fom_velo_t, rom_state_t>::value,"");
}

TEST(rom_galerkin_meta, validProjectorImplicitStepping)
{
  static_assert
    (pressio::rom::galerkin::constraints::projector_implicit_stepping
     <Projector, fom_velo_t, dmat_t, rom_state_t, rom_jac_t>::value,"");
}
