
# Components

The pressio C++ library is divided into several components:

| Name                							| Brief Description 													 	  | Links                                                							  | Corresponding header(s)    |
| ------------------                            | ---------------                                                             | ------------                                                                	  |  						|
| @m_span{m-text m-success}mpl@m_endspan        | metaprogramming functionalities                                             | [Documentation](md_pages_components_mpl.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/mpl)         	  | `<pressio/mpl.hpp>`  |
| @m_span{m-text m-success}utils@m_endspan      | common functionalities<br/>e.g., I/O helpers, static constants, etc         | [Documentation](md_pages_components_utils.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/utils)        	  | `<pressio/utils.hpp>` |
| @m_span{m-text m-success}type_traits@m_endspan| traits/detection classes   												  | [Documentation](md_pages_components_type_traits.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/type_traits)   	  | `<pressio/type_traits.hpp>`|
| @m_span{m-text m-success}expressions@m_endspan| expression classes for useful abstractions (span, diagonal, subspan, etc.)  | [Documentation](md_pages_components_expressions.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/expressions)        | `<pressio/expressions.hpp>`|
| @m_span{m-text m-success}ops@m_endspan        | shared-memory and distributed linear algebra kernels specializations        | [Documentation](md_pages_components_ops.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/ops)                | `<pressio/ops.hpp>`|
| @m_span{m-text m-success}qr@m_endspan         | QR factorization functionalities                                            | [Documentation](md_pages_components_qr.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/qr)                 | `<pressio/qr.hpp>`|
| @m_span{m-text m-success}solvers_linear@m_endspan    | linear solvers (wrappers around existing TPLs) 					  | [Documentation](md_pages_components_linsolvers.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/solvers_linear)            | `<pressio/solvers_linear.hpp>`|
| @m_span{m-text m-success}solvers_nonlinear@m_endspan | non-linear solvers <br> (e.g., Newton-Raphson, Gauss-Newton, Levenberg-Marquardt) | [Documentation](md_pages_components_nonlinsolvers.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/solvers_nonlinear) | `<pressio/solvers_nonlinear.hpp>`|
| @m_span{m-text m-success}ode@m_endspan        | explicit methods <br/>implict methods <br/> all   | <br/><br/>[Documentation](md_pages_components_ode.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/ode)   | `<pressio/ode_explicit.hpp>`<br/> `<pressio/ode_implicit.hpp>` <br/> `<pressio/ode.hpp>` |
| @m_span{m-text m-success}rom@m_endspan        | Galerkin<br/> LSPG<br/> WLS<br/> all      | <br/><br/><br/>[Documentation](md_pages_components_rom.html)<br/>[Code](https://github.com/Pressio/pressio/tree/main/include/rom) | `<pressio/rom_galerkin.hpp>` <br/> `<pressio/rom_lspg.hpp>` <br/> `<pressio/rom_wls.hpp>` <br/> `<pressio/rom.hpp>` |

The top-down order used above is informative of the dependency structure.
At the bottom of the stack we have the `rom` package which requires all the others.
This structure has several benefits.
* Maintability: `pressio` can be more easily developed and maintained since its components depend on one another through well-defined public interfaces,
and appropriate namespaces are used to organize classes.

* Selective usability: this modular framework allows users, if needed, to leverage invidual components.
For instance, if a user needs/wants just the QR methods, they can simply use that package,
and all the dependencies on the others are enabled automatically.


@m_class{m-block m-warning}

@par
When you use functionalities from a specific package, you should just include
the corresponding header and the dependencies (based on the explanation above) are included automatically.
For example, if you want to use Galerkin, you just need `#include <pressio/rom_galerkin.hpp>`
because all the needed packages are automatically included. There is not need to manually include all of them yourself.
In the future, we might refine further the granularity of the headers to allow a finer control.


<!-- @m_class{m-block m-warning} -->
<!-- @par One header to include them all -->
<!-- If you want to access *all* functionalities, you can use: -->
<!-- ```cpp -->
<!-- #include "pressio.hpp" -->
<!-- ``` -->
