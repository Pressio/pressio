
#ifndef PRESSIO_TEST_ROM_LSPG_STEADY_DEFAULT_CORRECT_FOMS_HPP_
#define PRESSIO_TEST_ROM_LSPG_STEADY_DEFAULT_CORRECT_FOMS_HPP_

#include <gtest/gtest.h>
#include "../../custom_data_types.hpp"

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#endif

struct TrivialFomContTimeEigen
{
  using scalar_type    = double;
  using state_type     = Eigen::VectorXd;
  using velocity_type  = state_type;
  int N_ = {};

  TrivialFomContTimeEigen(int N): N_(N){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.cols());
    return A;
  }

  void velocity(const state_type & u, scalar_type time, velocity_type & f) const
  {
    // EXPECT_TRUE(u.size()==f.size());
    // EXPECT_TRUE(u.size()==N_);
    // for (auto i=0; i<f.rows(); ++i){
    //  f(i) = u(i) + 1.;
    // }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     scalar_type time,
                     OperandType & A) const
  {
    // A = B;
    // A.array() += 1.;
  }
};

// struct TrivialFomSteadyCustomTypes
// {
//   using scalar_type       = double;
//   using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
//   using velocity_type     = state_type;
//   int N_ = {};

//   TrivialFomSteadyCustomTypes(int N): N_(N){}

//   velocity_type createVelocity() const{ return velocity_type(N_); }

//   template<class OperandType>
//   OperandType createApplyJacobianResult(const OperandType & B) const
//   {
//     OperandType A(N_, B.extent(1));
//     return A;
//   }

//   void velocity(const state_type & u, velocity_type & r) const
//   {
//     EXPECT_TRUE(u.extent(0)==r.extent(0));
//     EXPECT_TRUE(u.extent(0)==N_);

//     for (auto i=0; i<r.extent(0); ++i){
//      r(i) = u(i) + 1.;
//     }
//   }

//   template<class OperandType>
//   void applyJacobian(const state_type & state,
//                      const OperandType & B,
//                      OperandType & A) const
//   {
//     A = B;
//     for (std::size_t i=0; i< A.extent(0); ++i){
//       for (std::size_t j=0; j< A.extent(1); ++j){
//         A(i,j) += 1.;
//       }
//     }
//   }
// };

#endif
