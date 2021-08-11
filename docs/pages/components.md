
# Components

The pressio C++ library is divided into several components:

| Name                							| Brief Description 													 	  | Link                                                							  | Reference header(s)    |
| ------------------                            | ---------------                                                             | ------------                                                                	  |  						|
| @m_span{m-text m-success}mpl@m_endspan        | metaprogramming functionalities                                             | [Source](https://github.com/Pressio/pressio/tree/main/include/mpl)          	  | `#include<pressio_mpl.hpp>`  |
| @m_span{m-text m-success}utils@m_endspan      | common functionalities<br/>e.g., I/O helpers, static constants, etc         | [Source](https://github.com/Pressio/pressio/tree/main/include/utils)        	  | `#include<pressio_utils.hpp>` |
| @m_span{m-text m-success}type_traits@m_endspan| traits/detection classes   												  | [Source](https://github.com/Pressio/pressio/tree/main/include/type_traits)   	  | `#include<pressio_type_traits.hpp>`|
| @m_span{m-text m-success}expressions@m_endspan| expression classes for useful abstractions (span, diagonal, subspan, etc.)  | [Source](https://github.com/Pressio/pressio/tree/main/include/expressions)        | `#include<pressio_expressions.hpp>`|
| @m_span{m-text m-success}ops@m_endspan        | shared-memory and distributed linear algebra kernels specializations        | [Source](https://github.com/Pressio/pressio/tree/main/include/ops)                | `#include<pressio_ops.hpp>`| 
| @m_span{m-text m-success}qr@m_endspan         | QR factorization functionalities                                            | [Source](https://github.com/Pressio/pressio/tree/main/include/qr)                 | `#include<pressio_qr.hpp>`|
| @m_span{m-text m-success}solvers_linear@m_endspan    | linear solvers (wrappers around existing TPLs) 					  | [Source](https://github.com/Pressio/pressio/tree/main/include/solvers_linear)            | `#include<pressio_solvers_linear.hpp>`|
| @m_span{m-text m-success}solvers_nonlinear@m_endspan | non-linear solvers <br> (e.g., Newton-Raphson, Gauss-Newton, Levenberg-Marquardt) | [Source](https://github.com/Pressio/pressio/tree/main/include/solvers_nonlinear) | `#include<pressio_solvers_nonlinear.hpp>`|
| @m_span{m-text m-success}ode@m_endspan        | explicit only methods <br/>implict only methods <br/> all   | <br/><br/>[Source](https://github.com/Pressio/pressio/tree/main/include/ode)   | `#include<pressio_ode_explicit.hpp>`<br/> `#include<pressio_ode_implicit.hpp>` <br/> `#include<pressio_ode.hpp>` |
| @m_span{m-text m-success}rom@m_endspan        | Galerkin ROMs <br/> LSPG ROMs <br/> WLS ROMs <br/> all      | <br/><br/><br/>[Source](https://github.com/Pressio/pressio/tree/main/include/rom) | `#include<pressio_rom_galerkin.hpp>` <br/> `#include<pressio_rom_lspg.hpp>` <br/> `#include<pressio_rom_wls.hpp>` <br/> `#include<pressio_rom.hpp>` |

The top-down order used above is informative of the dependency structure.
For example, every package depends on `mpl`. The `ops` package depends only on `mpl`, `utils`, `containers`.
At the bottom of the stack we have the `rom` package which requires all the others.

This structure of the framework has several benefits.
* Maintability: `pressio` can be more easily developed and maintained since its components depend on one another through well-defined public interfaces,
and appropriate namespaces are used to organize classes.

* Selective usability: this modular framework allows users, if needed, to leverage invidual components.
For instance, if a user needs/wants just the QR methods, they can simply use that package,
and all the dependencies on the others are enabled automatically.


@m_class{m-block m-warning}

@par
When you use functionalities from a specific package, you should just include
the corresponding header and the dependencies (based on the explanation above) are included automatically.
For example, if you want to do Galerkin with explicit time integration,
you just do `#include <pressio_rom_galerkin.hpp>` because all the needed
packages are automatically included. There is not need to manually include all of them yourself.
In the future, we might refine further the granularity of the headers to allow a finer control.


@m_class{m-block m-warning}

@par One header to include them all
If you want to access *all* functionalities, you can use:
```cpp
#include "pressio.hpp"
```
