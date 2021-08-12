
# Linear solvers

Defined in header `<pressio/solvers_linear.hpp>`

Public namespace: `pressio::linearsolvers`


## Overview

First, we want to clearly state that this subpackage of pressio
does not provide linear solvers classes implemented from scratch
that work in general.
This is an intentional choice: there exist already many
open-source highly-optimized libraries that do so in both
a shared-memory, e.g. Armadillo, Eigen, Blaze, Magma, LAPACK,
and distributed environment, e.g. Trilinos, PETSc, HYPRE.

@m_span{m-text m-warning}The main goal of the pressio's linear solvers component
is to support usecases that are needed for ROMs.@m_endspan
As discussed in the introduction, ROMs involve "small"
*dense* systems that suitably fit on a *single compute node*.
Therefore, pressio's linear solvers are designed and implemented to be
wrappers (exposing a common interface) to existing shared-memory
linear solvers libraries.

Currently, we support solving linear systems that involve either
Eigen or Kokkos data structures. The reason for this is that Eigen
and Kokkos are the preferred choices (for now)
to implement ROMs operators. Kokkos, in particular, allows
to operate both on host and GPUs.
Later on, we can easily extend support to other libraries if needed.


## API
```cpp
pressio::linearsolvers::Solver<Tag, MatrixType>;
```

- `Tag`: tag type specifying which type of solver to use (more on this below)
- `MatrixType`: data type of the matrix of the linear system to solve
(used to choose the proper specialization).

We distinguish between iterative and direct solvers.
Obviously, depending on the structure of the matrix, some methods
might not be applicable or might not be efficient.

| Tag                         	                | Description   | Works for: 	|
|-----------------------------					|               |----------------	|
| `iterative::CG`               	|   Conjugate-gradient  | Eigen     	|
| `iterative::Bicgstab`         	|   Biconjugate gradient stabilized | Eigen     	|
| `iterative::LSCG`             	|   Conjugate-gradient  | Eigen     	|
| `direct::HouseholderQR`       	|   Uses on householder QR | Eigen     	|
| `direct::ColPivHouseholderQR` 	|   Uses on Householder QR with pivoting | Eigen     	|
| `direct::PartialPivLU`        	|   Uses LU factorization with partial pivoting | Eigen |
| `direct::potrsL`      |   Uses Cholesky, lower part	| Kokkos    |
| `direct::potrsU`      |   Uses Cholesky, upper part | Kokkos     	|
| `direct::getrs`       |   Uses LU factorization   | Kokkos     	|
| `direct::geqrf`       |   Uses QR fatorization  | Kokkos     	|


## Class Synopsis
```cpp
template<class MatrixType>
class LinearSolver
{
public:
  using matrix_type	= MatrixType;

  template<class VectorType>
  void solve(const MatrixType & A, const VectorType & b, VectorType & x);
};
```

## Example Usage

```cpp
//
// we want to solve `Ax = b`
//

using matrix_type = Eigen::MatrixXd;
using vector_type = Eigen::VectorXd;
matrix_type A(10, 10);
// fill A
vector_type b(10);
// fill b

namespace pls  = pressio::linearsolvers;
using tag	   = pls::direct::HouseholderQr;
using solver_t = pls::Solver<tag, matrix_type>;
solver mySolver;
vector_type x(10);
mySolver.solve(A, b, x);
```
