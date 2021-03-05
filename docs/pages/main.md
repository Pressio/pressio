# Pressio C++

*Leading-edge projection-based reduced order models (\proms) for
dynamical systems in science and engineering.*

You landed on the documentation of the C++ library,
which is one component of the [Pressio project](https://pressio.github.io/).


## In a nutshell

Pressio provides model reduction functionalities applicable to any system expressible as
@f[
\frac{d \boldsymbol{y}}{dt} =
\boldsymbol{f}(\boldsymbol{y},t; ...)
@f]
where @f$y@f$ is the full-order model (FOM) state,
@f$f@f$ is the FOM velocity, and @f$t@f$ is time.

We leverage this simple but expressive mathematical framework as a pivotal
design choice to enable a minimal application programming interface (API)
that is natural to dynamical systems.
To interface your application to Pressio, you need to create an *adapter class*
to expose the right-hand-side of your system and optionally the action of its Jacobian.
Schematically, this is shown below:
@image html frontpageschem.svg width=85%


## Getting Started

* learn about the [installation process](./md_pages_getstarted_build_and_install.html)

* read about the [packages](./md_pages_getstarted_packages.html) composing this C++ library

* explore the [tutorials](./md_pages_tutorials_tutorial1.html)

* dig into the [test subdirectory](https://github.com/Pressio/pressio/tree/master/tests/rom/burgers1d) of the C++ library


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
