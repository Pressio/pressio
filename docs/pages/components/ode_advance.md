
# ode advancers

Defined in header: `<pressio/ode_advancers.hpp>`

Public namespace: `pressio::ode`


## Overview

This component provides functionalities for "advancing" generic objects.
The main utility of these functions is for doing time integration
where one has a "stepper" and wants to use alternative strategies for advancing in time.

## API Synopsis:

```cpp
template<
  class StepperType,
  class StateType,
  class TimeType,
  class ...Args
  >
void advance_n_steps(StepperType & stepper,
					 StateType & state,
					 const TimeType start_time,
					 const TimeType time_step_size,
					 const ::pressio::ode::step_count_type num_steps,
					 Args && ... args);

template<
  class StepperType,
  class StateType,
  class TimeType,
  class StepSizeSetterType,
  class ...Args
  >
void advance_n_steps(StepperType & stepper,
					 StateType & state,
					 const TimeType start_time,
					 StepSizeSetterType && time_step_size_manager,
					 const ::pressio::ode::step_count_type num_steps,
					 Args && ... args);

template<
  class StepperType,
  class StateType,
  class TimeType,
  class ObserverType,
  class ...Args
  >
void advance_n_steps_and_observe(StepperType & stepper,
								 StateType & state,
								 const TimeType start_time,
								 const TimeType time_step_size,
								 const ::pressio::ode::step_count_type num_steps,
								 ObserverType & observer,
								 Args && ... args);

template<
  class StepperType,
  class StateType,
  class TimeType,
  class StepSizeSetterType,
  class ObserverType,
  class ...Args
  >
void advance_n_steps_and_observe(StepperType & stepper,
								 StateType & state,
								 const TimeType start_time,
								 StepSizeSetterType && time_step_size_manager,
								 const ::pressio::ode::step_count_type num_steps,
								 ObserverType & observer,
								 Args && ... args);
```

## Parameters

- `StepperType`:
  - represents the type of the object that can be stepped forward

- `StateType`:
  - type of the data structure you use for the state

- `TimeType`:
  - type of time: typically it is same as the value_type of the state

- `StepSizeSetterType`
  - type of object responsible for setting the time step size

- `ObserverType`:
  - type of the class responsible for "observing" the state during the
  time integration and potentially collect necessary data/metrics/statistics

- `Args...`:
  - any arguments needed to be passed to `doStep` method of the stepper


## Requirements

- The stepper class must conform to the following API:
  ```cpp
  class SteppableClass
  {
    public:
	void operator()(const StateType & current_state,
					const TimeType current_time,
					const TimeType time_step_size_to_use,
					const int32_t step_count
					/*[, other arguments needed: these will be forwarded to your
					stepper object by the advance methods as shown above ]*/);
  };
  ```

- The observer class must conform to the following API:
  ```cpp
  class ObserverClass
  {
    public:
    void operator()(const int32_t & step_count,
				    const TimeType current_time,
					const StateType & current_state) const;
  };
  ```

- `StepSizeSetterType` must conform to the following API:
  ```cpp
  class TimeStepManager
  {
    public:
    void operator()(const int32_t step_count,
				    const TimeType time,
					TimeType & time_step_size_to_set) const;
  };
  ```
