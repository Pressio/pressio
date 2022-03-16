
# Nonlinear Solvers - General Info


@m_class{m-note m-default}

@parblock
Defined in header `<pressio/solvers_nonlinear.hpp>`

Public namespace: `pressio::nonlinearsolvers`
@endparblock


## Overview

The pressio nonlinear solvers have been developed from the ground up using generic programming.
It is intended to be a library that can be used independently and with support for arbitrary data types.
It is not fully complete, and it can obviously be extended.
If you are interested in contributing or you need other methods, open an issue and/or PR.


Currently, we support the following algorithms:

| Name                                 | Doc                                                   | Purpose:                                                                                                                                                                                               |
|--------------------------------------|-------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Newton-Raphson                       | [Link](./md_pages_components_nonlinsolvers_nr.html)   | Systems of nonlinear equations (see e.g. [link](https://link.springer.com/content/pdf/bbm%3A978-3-319-69407-8%2F1.pdf), [link](https://www.cmu.edu/math/undergrad/suami/pdfs/2014_newton_method.pdf) ) |
| Gauss-Newton                         | [Link](./md_pages_components_nonlinsolvers_gn.html)   | Nonlinear least-squares problem.            (see [link](https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm) )                                                                                |
| Levenberg–Marquardt                  | [Link](./md_pages_components_nonlinsolvers_lm.html)   | Nonlinear least-squares problem.             (see [link](https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) )                                                                        |
| Iteratively reweighted least squares (irls) | [Link](./md_pages_components_nonlinsolvers_irls.html) | optimization problem formulated in a p-norm (see [link](https://en.wikipedia.org/wiki/Iteratively_reweighted_least_squares) )                                                                          |

<br/>

### A glimpse of the API

To create a solver, we expose specific factory functions for each algorithm.

```cpp
namespace pnonls = pressio::nonlinearsolvers;
auto solver = pnonls::create_newton_raphson(....);
auto solver = pnonls::create_gauss_newton(....);
```

Please refer to each method's documentation for the details on the API and arguments.

<br/>

## What are you responsible of?

If you want to use the pressio nonlinear solvers, you are responsible of:

- A: providing a problem object
- B: selecting (if needed) the convergence criterion
- C: selecting (if needed) the update method


### A: problem class

The problem class is what you use to *define your problem*.
An instance of a problem class is what you provide to pressio to
compute the needed operators to operate on. Obivously, this problem class must meet a specific API.

The pressio nonlinear solvers define four main types of problem APIs,
namely the "residual-jacobian", "fused residual-jacobian", "hessian-gradient" ,
and "fused hessian-gradient".
These problem concepts are discussed on [this page](md_pages_components_nonlinsolvers_system_api.html).
Depending on which concepts your problem meets, you can access certain algorithms but not
necessarily others.
The following table helps clarifying the problem API/algorithms admissibility:

| Algorithm                         | Residual-Jacobian API | Fused Residual-Jacobian API | Hessian-Gradient API | Fused Hessian-Gradient API |
|-----------------------------------|-----------------------|-----------------------------|----------------------|----------------------------|
| Newton-Raphson                    | admissible            | admissible                  | -                    | -                          |
| Gauss-Newton via Normal-Equations | admissible            | admissible                  | admissible           | admissible                 |
| Gauss-Newton via QR               | admissible            | admissible                  | -                    | -                          |
| Levenberg-Marquardt               | admissible            | admissible                  | admissible           | admissible                 |
| irls                              | admissible            | admissible                  | -                    | -                          |

Please refer to each method's documentation for the details on how you provide
a problem instance to pressio.



### B: Convergence and tolerance

The convergence criterion and associated tolerance are used to decide
why and when the solver needs to terminate.
We currently support these termination criteria:

| Enum value                                       | Description      | Currently supported for: |
|--------------------------------------------------|------------------|--------------------------|
| `Stop::AfterMaxIters`                            | self-explanatory | all algorithms           |
| `Stop::WhenCorrectionAbsoluteNormBelowTolerance` | self-explanatory | all algorithms           |
| `Stop::WhenCorrectionRelativeNormBelowTolerance` | self-explanatory | all algorithms           |
| `Stop::WhenResidualAbsoluteNormBelowTolerance`   | self-explanatory | all algorithms           |
| `Stop::WhenResidualRelativeNormBelowTolerance`   | self-explanatory | all algorithms           |
| `Stop::WhenGradientAbsoluteNormBelowTolerance`   | self-explanatory | least-squares solvers    |
| `Stop::WhenGradientRelativeNormBelowTolerance`   | self-explanatory | least-squares solvers    |


which you set/query using the following methods:

```cpp
class Solver{
  ...

  // set/query stopping criterion
  void setStoppingCriterion(Stop value);
  Stop stoppingCriterion() const;

  // set/query max number of iterations
  void setMaxIterations(iteration_type maxIters);
  iteration_type maxIterations() const;

  // this is used to set a single tol for all
  void setTolerance(scalar_type tolerance){ tolerances_.fill(std::move(tolerance)); }

  // finer-grained methods for setting tolerances
  void setCorrectionAbsoluteTolerance(scalar_type value);
  void setCorrectionRelativeTolerance(scalar_type value);
  void setResidualAbsoluteTolerance(scalar_type value);
  void setResidualRelativeTolerance(scalar_type value);
  void setGradientAbsoluteTolerance(scalar_type value);
  void setGradientRelativeTolerance(scalar_type value);

  // querying tolerances
  scalar_type correctionAbsoluteTolerance()const;
  scalar_type correctionRelativeTolerance()const;
  scalar_type residualAbsoluteTolerance()const;
  scalar_type residualRelativeTolerance()const;
  scalar_type gradientAbsoluteTolerance()const;
  scalar_type gradientRelativeTolerance()const;
  ...
};
```

### C: Setting the update

The update stage represents the *how* the current correction term is combined
with state to update the latter. We currently support the following:

| Name         | Enum value            | Description                         | Currently supported for: |
|--------------|-----------------------|-------------------------------------|--------------------------|
| Default      | `Update::Standard`    | @f$x_{n+1} = x_{n} + \lambda_{n}@f$ | all algorithms           |
| Armijo       | `Update::Armijo`      | [see this]()                        | Gauss-Newton             |
| LM-schedule1 | `Update::LMSchedule1` | [see this]()                        | Levenberg–Marquardt      |
| LM-schedule2 | `Update::LMSchedule2` | [see this]()                        | Levenberg–Marquardt      |

where @f$\lambda_{n}@f$ is the correction computed at the n-th iteration of the solver.

To set or query the update method, you use the following methods of the solver class:

```cpp
class Solver{
  ...
  void setUpdatingCriterion(Update value);
  Update updatingCriterion() const;
  ...
};
```


@m_class{m-block m-success}

@par By default, a nonlinear solver uses:
- update: `Update::Standard`;
- stopping: `Stop::WhenCorrectionAbsoluteNormBelowTolerance`;
- max number of iterations = 100
- tolerance = 0.000001 (for everything)


## A note on the solvers' design

If you are interested, here we provide a brief note describing the our design idea.
The design of the nonlinear solvers has been based on recognizing that, at a very high level,
a nonlinear solver operates by repeatedly updating a given "state" until a certain criterion is met,
and each "iteration" involves the following stages:

- A: computing/updating the operators

- B: computing the new correction term

- C: assessing convergence

- D: updating the state using the correction

This view forms the basis of our design approach: when a solver object is instantiated,
depending on the chosen algorithm, pressio behind the scenes instantiates the proper
classes needed for each of the stages above,
and properly composes them to instantiate the desired solver object.
