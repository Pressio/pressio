
# Pressio C++ Library


@m_class{m-frame m-default}

@parblock
*Advancing reduced order models (ROMs) for dynamical systems in science and engineering.*

This is the documentation of the [C++ library](https://github.com/Pressio/pressio), one component of the [Pressio ecosystem](https://pressio.github.io/).
@endparblock

<br/>

| Name                                                 | Description/Content                                                                        | Links                                                                                                                                                                                                                                                                                                                                                                       | Corresponding header(s)                                                                                                                                                            |
| ------------------                                   | ---------------                                                                            | ------------                                                                                                                                                                                                                                                                                                                                                                |                                                                                                                                                                                    |
| @m_span{m-text m-success}mpl@m_endspan               | metaprogramming functionalities                                                            | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio/mpl) <br/> [Documentation](md_pages_components_mpl.html)                                                                                                                                                                                                                                                        | `<pressio/mpl.hpp>`                                                                                                                                                                |
| @m_span{m-text m-success}utils@m_endspan             | logging, constants, various helpers, etc                                                | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio/utils)<br/>[Documentation](md_pages_components_utils.html)                                                                                                                                                                                                                                                      | `<pressio/utils.hpp>`                                                                                                                                                              |
| @m_span{m-text m-success}type_traits@m_endspan       | traits/detection classes                                                                   | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio/type_traits)<br/>[Documentation](md_pages_components_type_traits.html)                                                                                                                                                                                                                                          | `<pressio/type_traits.hpp>`                                                                                                                                                        |
| @m_span{m-text m-success}expressions@m_endspan       | classes for various abstractions (span, diagonal, subspan, etc.)                           | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio/expressions)<br/>[Documentation](md_pages_components_expressions.html)                                                                                                                                                                                                                                          | `<pressio/expressions.hpp>`                                                                                                                                                        |
| @m_span{m-text m-success}ops@m_endspan               | specializations of shared-memory and distributed linear algebra kernels                    | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio/ops)<br/>[Documentation](md_pages_components_ops.html)                                                                                                                                                                                                                                                          | `<pressio/ops.hpp>`                                                                                                                                                                |
| @m_span{m-text m-success}qr@m_endspan                | QR factorization functionalities                                                           | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio/qr)<br/>[Documentation](md_pages_components_qr.html)                                                                                                                                                                                                                                                            | `<pressio/qr.hpp>`                                                                                                                                                                 |
| @m_span{m-text m-success}solvers_linear@m_endspan    | linear solvers (wrappers around existing TPLs)                                             | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio/solvers_linear)<br/>[Documentation](md_pages_components_linsolvers.html)                                                                                                                                                                                                                                        | `<pressio/solvers_linear.hpp>`                                                                                                                                                     |
| @m_span{m-text m-success}solvers_nonlinear@m_endspan | <br/> general info <br/> Newton-Raphson <br/> Gauss-Newton <br/> Levenberg-Marquardt <br/> | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio/solvers_nonlinear) <br/> [Documentation](md_pages_components_nonlinsolvers_general.html) <br/> [Documentation](md_pages_components_nonlinsolvers_nr.html) <br/> [Documentation](md_pages_components_nonlinsolvers_gn.html) <br/> [Documentation](md_pages_components_nonlinsolvers_lm.html)                     | `<pressio/solvers_nonlinear.hpp>`                                                                                                                                                  |
| @m_span{m-text m-success}ode@m_endspan               | <br/> explicit steppers <br/>implicit steppers <br/> advancers <br/>                        | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio) <br/> [Documentation](md_pages_components_ode_steppers_explicit.html)<br/> [Documentation](md_pages_components_ode_steppers_implicit.html) <br/>[Documentation](md_pages_components_ode_advance.html)                                                                                                      | <br/> `<pressio/ode_steppers_explicit.hpp>` <br/> `<pressio/ode_steppers_implicit.hpp>`<br/> `<pressio/ode_advancers.hpp>` <br/> (for everything: `<pressio/ode.hpp>`)              |
| @m_span{m-text m-success}rom@m_endspan               | <br/>general info <br/> decoder <br/> Galerkin<br/> LSPG: steady<br/> LSPG: unsteady<br/> WLS<br/>                     | [Code](https://github.com/Pressio/pressio/tree/develop/include/pressio/rom) <br/>[Documentation](md_pages_components_rom_general.html) <br/>[Documentation](md_pages_components_rom_decoder.html) <br/> [Documentation](md_pages_components_rom_galerkin.html) <br/> [Documentation](md_pages_components_rom_lspg_steady.html) <br/> [Documentation](md_pages_components_rom_lspg_unsteady.html) <br/>  [Documentation](md_pages_components_rom_wls.html) <br/> | <br/> <br/> `<pressio/rom_decoder.hpp>` <br/> `<pressio/rom_galerkin.hpp>` <br/> `<pressio/rom_lspg.hpp>` <br/> `<pressio/rom_lspg.hpp>`<br/> `<pressio/rom_wls.hpp>` <br/> (for everything: `<pressio/rom.hpp>`) |

<!-- This structure benefits: -->
<!-- * Maintability: the components depend on one another through well-defined public interfaces, -->
<!-- and appropriate namespaces are used to properly scope the code. -->

<!-- * Selective usability: users can leverage invidual functionalities. -->
<!-- This allows finer-grained control on what you include and use. -->

## Get Started

* [read the introduction](./md_pages_introduction.html) providing an overview, objectives and design ideas

* [how to install](./md_pages_installation.html): it is a header-only library, should be trivial

* [explore the tutorials](https://pressio.github.io/pressio-tutorials/html/index.html)


<!-- ## What if your types are not natively supported in pressio? -->

<!-- Check if your types are supported by lookig at the -->
<!-- [dependencies](md_pages_getstarted_build_and_install.html): if they are -->
<!-- listed there, most likely you are good to go, and you don't need to provide extra information to pressio. -->

<!-- Not supported? You can file an [issue](https://github.com/Pressio/pressio/issues) -->
<!-- to request it and wait on it, or can proceed -->
<!-- as in [tutorialsB](./md_pages_tutorials_tutorial1udops.html). Or do both! -->


## License and Citation
The full license (BSD-3) is available [here](https://pressio.github.io/various/license/).

We are working on publishing this: you can find our arXiv preprint at: https://arxiv.org/abs/2003.07798

## Questions?
Find us on Slack: https://pressioteam.slack.com or
open an issue on [github](https://github.com/Pressio/pressio).


<!--
@m_class{m-note m-success}

Pressio is an open-source project aimed at enabling leading-edge projection-based
reduced order models (\proms) for dynamical systems in science and engineering.

## Motivation
Projection-based model reduction refers to a class of surrogate models
that reduce the number of degrees
of freedom in the full-order model (FOM) through a projection process.
This projection step applied to the governing equations often enables one
to make stronger performance guarantees
(e.g., of structure preservation, of accuracy via adaptivity) than other
surrogates like data-fits and perform more accurate *a posteriori*
error analysis (e.g., via *a posteriori* error bounds or error models).

Despite these benefits, the practical challenges of
implementing model-reduction techniques in large-scale codes often
precludes their adoption in practice; this occurs because standard implementations
require modifying low-level operations and solvers for each simulation code of interest.
This implementation strategy is not practical or sustainable
in many modern settings, because industrial simulation codes often evolve rapidly,
institutions may employ dozens of simulation codes for different analyses,
and commercial codes typically do not expose the required low-level
operators and solvers.


@m_class{m-note m-success}

Pressio aims to mitigate the implementation burden of projection-based model
reduction in large-scale applications without compromising performance.


## Main steps of pROMs
Projection-based model reduction can be broken into three main steps,
namely data collection, basis creation, and ROM deployment.

- data collection: \todo (all)

- compute basis: \todo (all)

- create/run the ROM: \todo (all)


@m_class{m-block m-warning}

@par
pressioproj currently contains capabilities to perform the last step.
\todo Say that we have plans for the other steps too.
Maybe at some point we will provide tools to run the samples,
but for now that is not a huge priority. we can develop something
later on to aid this step. For example interfacing with efficient
POD libraries, providing tools for specific mesh formats (exodus).
 -->

<!--
## The Pressio framework
\pressioproj is a computational *framework*, comprising a (growing) collection of repositories :

* [pressio](https://github.com/Pressio/pressio): &emsp;&ensp;&emsp;&emsp;&ensp;core C++ library based on generic programming;

<!-- to support applications with arbitrary data types; -->
<!-- [pressio4py](https://github.com/Pressio/pressio4py): &emsp;&emsp;&nbsp;&nbsp;Python bindings for the core Pressio C++ functionalities; -->
<!-- [pressio-builder](https://github.com/Pressio/pressio-builder): &nbsp;&nbsp;&nbsp;auxiliary bash scripts for building/testing; -->
<!-- [pressio-tutorials](https://github.com/Pressio/pressio-tutorials): &nbsp;tutorials explaining how to use `pressio` and its functionalities.

## Where to go from here
If you are new and want to learn more, start from the [userguide](./md_pages_get_started.html)
and see how to install and use pressio, or you can jump directly
to the [tutorials](./md_pages_tutorials.html)
and/or [examples](md_pages_examples.html) -->
