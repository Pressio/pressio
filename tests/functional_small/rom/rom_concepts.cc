
#include <gtest/gtest.h>
#include "pressio/rom_decoder.hpp"

#include "pressio/rom/predicates/all.hpp"
#include "pressio/rom/constraints/all.hpp"

namespace
{
struct ValidFomState{
  ValidFomState(const ValidFomState &) = default;
};

struct InvalidFomState{
  InvalidFomState(const InvalidFomState &) = delete;
};

template <class MatrixType, class FomStateType>
struct ValidDecoder
{
  using jacobian_type  = MatrixType;
  using fom_state_type = FomStateType;

public:
  template <typename RomStateType>
  void applyMapping(const RomStateType & romState, fom_state_type & result) const;
  template <typename RomStateType> void updateJacobian(const RomStateType &);
  const jacobian_type & jacobianCRef() const;
};

struct MyArbFomVector{};
struct MyArbFomMatrix{};

struct ValidSteadyFomSystem{
  using scalar_type = double;
  using state_type  = MyArbFomVector;
  using residual_type = state_type;
public:
  residual_type createResidual() const;
  MyArbFomMatrix createApplyJacobianResult(const MyArbFomMatrix &) const;
  void residual(const state_type &, residual_type &) const;
  void applyJacobian(const state_type &, const MyArbFomMatrix &, MyArbFomMatrix &) const;
};

struct ValidFomSystemOnlyVelocity{
  using scalar_type = double;
  using state_type  = MyArbFomVector;
  using velocity_type = state_type;
public:
  velocity_type createVelocity() const;
  void velocity(const state_type &,const scalar_type &,velocity_type &) const;
};

struct ValidFomSystemVelocityAndApplyJacobian{
  using scalar_type = double;
  using state_type  = MyArbFomVector;
  using velocity_type = state_type;
public:
  velocity_type createVelocity() const;
  MyArbFomMatrix createApplyJacobianResult(const MyArbFomMatrix & B) const;

  void velocity(const state_type &,const scalar_type &,velocity_type &) const;
  void applyJacobian(const state_type &,const MyArbFomMatrix &,const scalar_type &, MyArbFomMatrix &) const;
};

struct ValidFomSystemDiscreteTimeTwoStates{
  using scalar_type                 = double;
  using state_type                  = MyArbFomVector;
  using discrete_time_residual_type = state_type;

public:
  discrete_time_residual_type createDiscreteTimeResidual() const;
  MyArbFomMatrix createApplyDiscreteTimeJacobianResult(const MyArbFomMatrix & B) const;

  template <typename step_t>
  void discreteTimeResidual(const step_t & step,
            const scalar_type & time,
            const scalar_type & dt,
            discrete_time_residual_type & R,
            const state_type &,
            const state_type &) const;

  template <typename step_t>
  void applyDiscreteTimeJacobian(const step_t & step,
           const scalar_type & time,
           const scalar_type & dt,
           const MyArbFomMatrix & B,
           MyArbFomMatrix & A,
           const state_type &,
           const state_type &) const;
};

struct ValidFomSystemDiscreteTimeThreeStates{
  using scalar_type                 = double;
  using state_type                  = MyArbFomVector;
  using discrete_time_residual_type = state_type;

public:
  discrete_time_residual_type createDiscreteTimeResidual() const;
  MyArbFomMatrix createApplyDiscreteTimeJacobianResult(const MyArbFomMatrix & B) const;

  template <typename step_t>
  void discreteTimeResidual(const step_t & step,
            const scalar_type & time,
            const scalar_type & dt,
            discrete_time_residual_type & R,
            const state_type &,
            const state_type &,
            const state_type &) const;

  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
           const scalar_type & time,
           const scalar_type & dt,
           const MyArbFomMatrix & B,
           MyArbFomMatrix & A,
           const state_type &,
           const state_type &,
           const state_type &) const;
};

struct ValidFomSystemDiscreteTimeFourStates{
  using scalar_type                 = double;
  using state_type                  = MyArbFomVector;
  using discrete_time_residual_type = state_type;

public:
  discrete_time_residual_type createDiscreteTimeResidual() const;
  MyArbFomMatrix createApplyDiscreteTimeJacobianResult(const MyArbFomMatrix & B) const;

  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
            const scalar_type & time,
            const scalar_type & dt,
            discrete_time_residual_type & R,
            const state_type &,
            const state_type &,
            const state_type &,
            const state_type &) const;

  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
           const scalar_type & time,
           const scalar_type & dt,
           const MyArbFomMatrix & B,
           MyArbFomMatrix & A,
           const state_type &,
           const state_type &,
           const state_type &,
           const state_type &) const;
};

}//end namespace anonim

TEST(rom, concepts_fom_state)
{
  namespace prom = pressio::rom;
  static_assert(prom::fom_state<ValidFomState>::value, "");
  static_assert(!prom::fom_state<InvalidFomState>::value, "");
}

TEST(rom, concepts_decoder)
{
  namespace prom = pressio::rom;
  using decoder_t = ValidDecoder<MyArbFomMatrix, MyArbFomVector>;
  using rom_state_t = std::vector<double>;
  static_assert(prom::decoder<decoder_t, rom_state_t>::value, "");
}

TEST(rom, concepts_fom_system_steady)
{
  namespace prom = pressio::rom;

  ///////////////////////////////////
  // a steady FOM class should
  {
    using fom_t = ValidSteadyFomSystem;

    // - satisfy the steady concept
    static_assert(prom::steady_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert(prom::most_likely_steady_fom_system<
      fom_t>::value, "");

    // - NOT satisfy all the others
    static_assert( !prom::most_likely_continuous_time_fom_system<
      fom_t>::value, "");
    static_assert( !prom::most_likely_discrete_time_fom_system<
      fom_t>::value, "");
    static_assert( !prom::continuous_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::continuous_time_fom_system_without_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::continuous_time_fom_system_with_at_least_velocity<
      fom_t>::value, "");
    static_assert( !prom::continuous_time_fom_system<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::discrete_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, 1, MyArbFomMatrix>::value, "");
    static_assert( !prom::discrete_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, 2, MyArbFomMatrix>::value, "");
  }

  ///////////////////////////////////
  // a FOM class with only velocty should
  {
    using fom_t = ValidFomSystemOnlyVelocity;

    // - satisfy these
    static_assert(prom::most_likely_continuous_time_fom_system<
      fom_t>::value, "");
    static_assert(prom::continuous_time_fom_system_with_at_least_velocity<
      fom_t>::value, "");
    static_assert(prom::continuous_time_fom_system_without_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert(prom::continuous_time_fom_system<
      fom_t, MyArbFomMatrix>::value, "");

    // fail all these
    static_assert( !prom::continuous_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::steady_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::most_likely_steady_fom_system<
      fom_t>::value, "");
    static_assert( !prom::most_likely_discrete_time_fom_system<
      fom_t>::value, "");
    static_assert( !prom::discrete_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, 1, MyArbFomMatrix>::value, "");
    static_assert( !prom::discrete_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, 2, MyArbFomMatrix>::value, "");
  }

  ////////////////////////////////////////////////
  // a FOM class with velocty and apply jac should
  {
    using fom_t = ValidFomSystemVelocityAndApplyJacobian;

    // - satisfy these
    static_assert(prom::most_likely_continuous_time_fom_system<
      fom_t>::value, "");
    static_assert(prom::continuous_time_fom_system_with_at_least_velocity<
      fom_t>::value, "");
    static_assert(prom::continuous_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert(prom::continuous_time_fom_system<
      fom_t, MyArbFomMatrix>::value, "");

    // fail all these
    static_assert( !prom::continuous_time_fom_system_without_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::steady_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::most_likely_steady_fom_system<
      fom_t>::value, "");
    static_assert( !prom::most_likely_discrete_time_fom_system<
      fom_t>::value, "");
    static_assert( !prom::discrete_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, 1, MyArbFomMatrix>::value, "");
    static_assert( !prom::discrete_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, 2, MyArbFomMatrix>::value, "");
  }

  ////////////////////////////////////////////////
  // a FOM class with discrete api 2 states should
  {
    using fom_t = ValidFomSystemDiscreteTimeTwoStates;

    // - satisfy these
    static_assert(prom::most_likely_discrete_time_fom_system<
      fom_t>::value, "");
    static_assert(prom::discrete_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, 2, MyArbFomMatrix>::value, "");

    // fail all these
    static_assert( !prom::most_likely_continuous_time_fom_system<
      fom_t>::value, "");
    static_assert( !prom::continuous_time_fom_system_with_at_least_velocity<
      fom_t>::value, "");
    static_assert( !prom::continuous_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::continuous_time_fom_system<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::continuous_time_fom_system_without_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::steady_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::most_likely_steady_fom_system<
      fom_t>::value, "");
  }


  ////////////////////////////////////////////////
  // a FOM class with discrete api 3 states should
  {
    using fom_t = ValidFomSystemDiscreteTimeThreeStates;

    // - satisfy these
    static_assert(prom::most_likely_discrete_time_fom_system<
      fom_t>::value, "");
    static_assert(prom::discrete_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, 3, MyArbFomMatrix>::value, "");

    // fail all these
    static_assert( !prom::most_likely_continuous_time_fom_system<
      fom_t>::value, "");
    static_assert( !prom::continuous_time_fom_system_with_at_least_velocity<
      fom_t>::value, "");
    static_assert( !prom::continuous_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::continuous_time_fom_system<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::continuous_time_fom_system_without_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::steady_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::most_likely_steady_fom_system<
      fom_t>::value, "");
  }

  ////////////////////////////////////////////////
  // a FOM class with discrete api 4 states should
  {
    using fom_t = ValidFomSystemDiscreteTimeFourStates;

    // - satisfy these
    static_assert(prom::most_likely_discrete_time_fom_system<
      fom_t>::value, "");
    static_assert(prom::discrete_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, 4, MyArbFomMatrix>::value, "");

    // fail all these
    static_assert( !prom::most_likely_continuous_time_fom_system<
      fom_t>::value, "");
    static_assert( !prom::continuous_time_fom_system_with_at_least_velocity<
      fom_t>::value, "");
    static_assert( !prom::continuous_time_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::continuous_time_fom_system<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::continuous_time_fom_system_without_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::steady_fom_system_with_user_provided_apply_jacobian<
      fom_t, MyArbFomMatrix>::value, "");
    static_assert( !prom::most_likely_steady_fom_system<
      fom_t>::value, "");
  }
}
