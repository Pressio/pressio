
# FOM Adapter API

<!-- As discussed in [this](./md_pages_getstarted_pressio_app.html) section of the get started page, -->
The adapter class sits between pressio and the external application.
This adapter class---if needed---must meet a specific API.

Pressio supports two variants of the API,
a so-called *continuous-time* and *discrete-time* version.
The continuous-time API operates such that the user is responsible
to compute the continuous-time operators, e.g., the velocity, and pressio assembles the
discrete-time operators. It is an API that very expressive of the formulation
on which pressio relies on.
The discrete-time API is designed such that the user is given the
necessary operators and operands needed to assemble the
time-discrete operators directly.

\todo: why are we forced to keep these two separate? explain

\todo: discuss pros and cons of the two APIs

The following table illustrates what API is supported for each ROM

| ROM Method                        | Continuous-time API support | Discrete-time API support |
| ------------------                | ---------------             | ---------------           |
| Galerkin (explicit time stepping) | yes                         | no                        |
| Galerkin (implicit time stepping) | yes                         | yes                       |
| LSPG                              | yes                         | yes                       |
| WLS                               | yes                         | yes                       |
