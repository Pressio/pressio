
# What is Pressio?

@m_class{m-block m-default} @par
@parblock
\pressioproj is an open-source project aimed at enabling
leading-edge projection-based reduced order models (\proms)
for large-scale nonlinear dynamical systems in science and engineering.
@endparblock

## Why Pressio?
Projection-based model reduction refers to a class of surrogate models
that reduce the number of degrees
of freedom in the high-fidelity model through a projection process.
This projection step applied to governing equations often enables one
to make stronger performance guarantees
(e.g., of structure preservation, of accuracy via adaptivity) than other
surrogates like data-fits and perform more accurate *a posteriori*
error analysis (e.g., via *a posteriori* error bounds or error models).

Despite these benefits, the practical challenges of
implementing nonlinear model-reduction techniques in large-scale codes often
precludes their adoption in practice; this occurs because standard implementations
require modifying low-level operations and solvers for each simulation code of interest.
This implementation strategy is simply not practical or sustainable
in many modern settings, because industrial simulation codes often evolve rapidly,
institutions may employ dozens of simulation codes for different analyses,
and commercial codes typically do not expose the required low-level
operators and solvers.

@m_span{m-text m-success}
\pressioproj aims to break this barrier by mitigating the implementation burden of
nonlinear model reduction in large-scale applications without compromising performance.
@m_endspan

## What steps of pROMs does Pressio cover?
Projection-based model reduction can be broken into three main steps,
namely data collection, basis creation, and ROM deployment.
\pressioproj currently contains capabilities to perform the last step.

- data collection: \todo (users do this. Maybe at some point we will provide tools
to run the samples, but for now that is not a huge priority.)

- compute basis: \todo (this is done externally. But we can develop something
later on to aid this step. For example interfacing with efficient
POD libraries, providing tools for specific mesh formats (exodus).)

- create/run the ROM: \todo


## What does Pressio include?
\pressioproj is a computational framework, comprising a (growing) collection of repositories :

* [pressio](https://github.com/Pressio/pressio): &emsp;&ensp;&emsp;&emsp;&ensp;core C++ library based on generic programming;
<!-- to support applications with arbitrary data types; -->
* [pressio4py](https://github.com/Pressio/pressio4py): &emsp;&emsp;&nbsp;&nbsp;Python bindings for the core Pressio C++ functionalities;
* [pressio-builder](https://github.com/Pressio/pressio-builder): &nbsp;&nbsp;&nbsp;auxiliary repo with bash scripts for building/testing.
* [pressio-tutorials](https://github.com/Pressio/pressio-tutorials): &nbsp;tutorials explaining how to use `pressio` and its functionalities;
