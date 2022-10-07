
#include <gtest/gtest.h>
#include "pressio/ode_concepts.hpp"

namespace{
  using namespace pressio::ode;

  struct Stepper1{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    void operator()(state_type &	    /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount		    /*unused*/,
		    StepSize<independent_variable_type> /*unused*/){}
  };

  struct Stepper2{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // wrong order
    void operator()(state_type &	    /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    StepCount		    /*unused*/){}
  };

  struct Stepper3{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // missing ref here
    void operator()(state_type              /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/){}
  };

  struct Stepper4{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // rvalue ref here
    void operator()(state_type &&           /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/){}
  };


  struct AuxThing1{};
  struct AuxThing2{};

  struct VarStepper1{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    void operator()(state_type &            /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/){}
  };

  struct VarStepper2{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // wrong order
    void operator()(state_type &            /*unused*/,
		    StepCount               /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/){}
  };

  struct VarStepper3{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // missing ref for state
    void operator()(state_type              /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/){}
  };

  struct VarStepper4{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // rvalue ref for state
    void operator()(state_type &&           /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/){}
  };


  struct VarStepper5{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    void operator()(state_type &            /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/,
		    AuxThing2               /*unused*/)
    {}
  };

} //end anonym namespace

TEST(ode, concepts__stepper)
{
  using namespace pressio::ode;

#ifdef PRESSIO_ENABLE_CXX20
  static_assert( Steppable<Stepper1>, "");
  static_assert(!Steppable<Stepper2>, "");
  static_assert(!Steppable<Stepper3>, "");
  static_assert(!Steppable<Stepper4>, "");
#else
  static_assert( Steppable<Stepper1>::value, "");
  static_assert(!Steppable<Stepper2>::value, "");
  static_assert(!Steppable<Stepper3>::value, "");
  static_assert(!Steppable<Stepper4>::value, "");
#endif
}

TEST(ode, concepts_variadic_stepper)
{
  using namespace pressio::ode;

#ifdef PRESSIO_ENABLE_CXX20
  static_assert(!Steppable<VarStepper1>, "");
  static_assert(!Steppable<VarStepper2>, "");
  static_assert(!Steppable<VarStepper3>, "");
  static_assert(!Steppable<VarStepper4>, "");
  static_assert(!Steppable<VarStepper5>, "");

  static_assert(SteppableWithAuxiliaryArgs<
		VarStepper1, AuxThing1>, "");
  static_assert(!SteppableWithAuxiliaryArgs<
		VarStepper2, AuxThing1>, "");
  static_assert(!SteppableWithAuxiliaryArgs<
		VarStepper3, AuxThing1>, "");
  static_assert(!SteppableWithAuxiliaryArgs<
		VarStepper4, AuxThing1>, "");

  static_assert(!SteppableWithAuxiliaryArgs<
		VarStepper5, AuxThing1>, "");
  static_assert(SteppableWithAuxiliaryArgs<
		VarStepper5, AuxThing1, AuxThing2>, "");
#else
  static_assert(!Steppable<VarStepper1>::value, "");
  static_assert(!Steppable<VarStepper2>::value, "");
  static_assert(!Steppable<VarStepper3>::value, "");
  static_assert(!Steppable<VarStepper4>::value, "");
  static_assert(!Steppable<VarStepper5>::value, "");

  static_assert(SteppableWithAuxiliaryArgs<
		void, VarStepper1, AuxThing1>::value, "");
  static_assert(!SteppableWithAuxiliaryArgs<
		void, VarStepper2, AuxThing1>::value, "");
  static_assert(!SteppableWithAuxiliaryArgs<
		void, VarStepper3, AuxThing1>::value, "");
  static_assert(!SteppableWithAuxiliaryArgs<
		void, VarStepper4, AuxThing1>::value, "");
  static_assert(!SteppableWithAuxiliaryArgs<
		void, VarStepper5, AuxThing1>::value, "");
  static_assert(SteppableWithAuxiliaryArgs<
		void, VarStepper5, AuxThing1, AuxThing2>::value, "");
#endif
}
