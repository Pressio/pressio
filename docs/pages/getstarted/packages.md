
# Packages

@m_class{m-block m-info}

@par What is this page about?
This page describes the structure of C++ pressio library.
By the end, you should be able to understand the structure
with its main packages, and the logic behind.


The pressio C++ library is divided into several packages:

| Package <br> name                             | Description                                                                                     | Link                                                                             |
| ------------------                            | ---------------                                                                                 | ------------                                                                     |
| @m_span{m-text m-success}mpl@m_endspan        | metaprogramming functionalities                                                                 | [Source](https://github.com/Pressio/pressio/tree/master/packages/mpl/src)        |
| @m_span{m-text m-success}utils@m_endspan      | common functionalities, e.g., I/O helpers, static constants, etc                                | [Source](https://github.com/Pressio/pressio/tree/master/packages/utils/src)      |
| @m_span{m-text m-success}containers@m_endspan | wrappers for vectors, matrices and multi-vectors, <br> expressions (span, diagonal and subspan) | [Source](https://github.com/Pressio/pressio/tree/master/packages/containers/src) |
| @m_span{m-text m-success}ops@m_endspan        | shared-memory and distributed linear algebra kernels                                            | [Source](https://github.com/Pressio/pressio/tree/master/packages/ops/src)        |
| @m_span{m-text m-success}apps@m_endspan       | suites of mini-apps used for basic testing                                                      | [Source](https://github.com/Pressio/pressio/tree/master/packages/apps/src)       |
| @m_span{m-text m-success}qr@m_endspan         | QR factorization functionalities                                                                | [Source](https://github.com/Pressio/pressio/tree/master/packages/qr/src)         |
| @m_span{m-text m-success}solvers@m_endspan    | linear and non-linear solvers <br> (e.g., Newton-Raphson, Gauss-Newton, Levenberg-Marquardt)    | [Source](https://github.com/Pressio/pressio/tree/master/packages/solvers/src)    |
| @m_span{m-text m-success}ode@m_endspan        | explicit and implict time steppers and integrators                                              | [Source](https://github.com/Pressio/pressio/tree/master/packages/ode/src)        |
| @m_span{m-text m-success}rom@m_endspan        | reduced-order modeling algorithms                                                               | [Source](https://github.com/Pressio/pressio/tree/master/packages/rom/src)        |

The top-down order used above is informative of the dependency structure.
For example, every package depends on `mpl`. The `ops` package depends only on `mpl`, `utils`, `containers`.
At the bottom of the stack we have the `rom` package which requires all the others.

Splitting the framework into separate packages has several benefits.
* Maintability: `pressio` can be more easily developed and maintained since packages depend on one another through well-defined public interfaces,
and appropriate namespaces are used to organize classes.

* Selective usability: this modular framework allows users, if needed, to leverage invidual packages.
For instance, if a user needs/wants just the QR methods, they can simply use that package,
and all the dependencies on the others are enabled automatically.
