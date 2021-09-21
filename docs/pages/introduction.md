
# Introduction

Finish: graphic needs to be updated, content to add.

- ROMs are dense
- pressio was designed with that in mind


## In a nutshell

Pressio can be applied to any dynamical system expressible in
a *continuous-time* form as
@f[
\frac{d \boldsymbol{y}}{dt} =
\boldsymbol{f}(\boldsymbol{y},t; ...)
@f]
and/or in a *discrete-time* form
@f[
\boldsymbol{R}(\boldsymbol{y}, \boldsymbol{y_{n-1}}, ..., t_n, dt_n; ...) = \boldsymbol{0}
@f]

Here, @f$y@f$ is the full-order model (FOM) state,
@f$f@f$ the FOM velocity, @f$t@f$ is time, and @f$R@f$ is the residual.

We leverage this expressive mathematical framework as a pivotal
design choice to enable a minimal application programming interface (API)
that is natural to dynamical systems: you choose the formulation
more convenient to you, and interface your application to
Pressio by creating a corresponding *adapter class* to expose
the operators needed for the chosen formulation.
In general, you don't need to support both: there are advantages and disadvantages for both,
and sometimes the choice is dictated directly by your native application (for example,
in some cases it might be easier to directly expose the residual).
Read [the doc page](md_pages_components_rom_fom_apis.html)
to learn more about the adapter classes and see code templates.

@image html frontpageschem.svg width=70%
