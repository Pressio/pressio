
\page frontpage Overview


@m_class{m-block m-success}

@par What is Pressio?
	An open-source project aimed at enabling leading-edge projection-based
	reduced order models (\proms) for dynamical systems in science and engineering.


## Motivation
Projection-based model reduction refers to a class of surrogate models
that reduce the number of degrees
of freedom in the high-fidelity model through a projection process.
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


@m_class{m-block m-success}

@par Objective
	To mitigate the implementation burden of projection-based model
	reduction in large-scale applications without compromising performance.


## Main steps of pROMs
Projection-based model reduction can be broken into three main steps,
namely data collection, basis creation, and ROM deployment.

- data collection: \todo (all)

- compute basis: \todo (all)

- create/run the ROM: \todo (all)

## What steps does Pressio currently cover?
\pressioproj currently contains capabilities to perform the last step.
\todo Say that we have plans for the other steps too.
Maybe at some point we will provide tools to run the samples,
but for now that is not a huge priority. we can develop something
later on to aid this step. For example interfacing with efficient
POD libraries, providing tools for specific mesh formats (exodus).


## The Pressio framework
\pressioproj is a computational *framework*, comprising a (growing) collection of repositories :

* [pressio](https://github.com/Pressio/pressio): &emsp;&ensp;&emsp;&emsp;&ensp;core C++ library based on generic programming;
<!-- to support applications with arbitrary data types; -->
* [pressio4py](https://github.com/Pressio/pressio4py): &emsp;&emsp;&nbsp;&nbsp;Python bindings for the core Pressio C++ functionalities;
* [pressio-builder](https://github.com/Pressio/pressio-builder): &nbsp;&nbsp;&nbsp;auxiliary bash scripts for building/testing;
* [pressio-tutorials](https://github.com/Pressio/pressio-tutorials): &nbsp;tutorials explaining how to use `pressio`
and its functionalities.

## Where to go from here
If you are new and want to learn more, start from the \ref userguide and see how to install and use pressio,
or you can jump directly to the \ref tutorials and/or \ref examples.
