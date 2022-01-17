
# ops

<!--
Enable TOC on sidebar for easy navigation ?

@tableofcontents
-->

@m_class{m-note m-default}

@parblock
Defined in header: `<pressio/ops.hpp>`

Public namespace: `pressio::ops`
@endparblock

## Overview

Operations provide algebraic routines and copy/assignment utilities that can be executed on containers and aim to offer relevant type and rank flexibility.

## Math

#### abs()

@m_class{m-block m-default}

@parblock
```cpp
template <class T1, class T2>
void abs(T1& y, const T2& x)
```
- Computes the absolute value of `x` storing the result in `y`: @f$y_i := |x_i|@f$
- Works on vectors (matrix support in progress).
@endparblock

#### min() and max()

@m_class{m-block m-default}

@parblock
```cpp
template <typename T>
Traits<T>::scalar_type min(const T& obj)

template <typename T>
Traits<T>::scalar_type max(const T& obj)
```
- `min()` selects the minimal and `max()` the maximal value from the elements contained in `obj`.
@endparblock

#### scale()

@m_class{m-block m-default}

@parblock
```cpp
template <typename T>
void scale(T& o, const Traits<T>::scalar_type value)
```
- Scales each element in `o` (vector or matrix) by scalar `value`: @f$o_\lambda := o_\lambda value@f$<br> (where @f$\lambda@f$ is rank one or rank two index). 
@endparblock

#### pow() and abs_pow()

@m_class{m-block m-default}

@parblock
```cpp
template <typename T>
void pow(T& x, const Traits<T>::scalar_type& exponent)
```
- Replaces each element in `x` with computed  expotential: @f$x_i := x_i^{expotent}@f$
@endparblock

@m_class{m-block m-default}

@parblock
```cpp
template <typename T1, typename T2>
void pow(T1& y, const T2& x, const Traits<T1>::scalar_type& exponent)
```
- Computes the expotential of `x`: @f$y_i := x_i^{expotent}@f$
@endparblock

@m_class{m-block m-default}

@parblock
```cpp
template <typename T1, typename T2>
void abs_pow(T1& y, const T2& x, const Traits<T>::scalar_type& exponent)
```
- Computes the expotential of absolute `x`: @f$y_i := |x_i|^{expotent}@f$
@endparblock

@m_class{m-block m-default}

@parblock
```cpp
template <typename T1, typename T2>
void abs_pow(
  T1& y,
  const T2& x,
  const Traits<T1>::scalar_type& exponent,
  const Traits<T1>::scalar_type& eps)
```
- Computes the *negative* expotential (requires `expotent < 0`) of absolute `x`: @f$y_i := \frac{1}{ max(|x_i|^{-expotent}, \epsilon) }@f$
@endparblock

#### norm1() and norm2()

@m_class{m-block m-default}

@parblock
```cpp
template <typename T>
Traits<T>::scalar_type norm1(const T& o)

template <typename T>
Traits<T>::scalar_type norm2(const T& o)
```
- `norm1()` computes and returns L1 norm for `o`: @f$||o||_1 = \sum_i{ |o_i| }@f$
- `norm2()` computes and returns L2 norm for `o`: @f$||o||_2 = \sqrt{ \sum_i{ o_i^2 }}@f$
@endparblock

#### add_to_diagonal()

@m_class{m-block m-default}

@parblock
```cpp
template <typename T>
void add_to_diagonal(T& mtx, Traits<T>::scalar_type value)
```
- Adds scalar `value` to each diagonal element of `mtx` (`mtx(i, i)`).
@endparblock

#### dot()

@m_class{m-block m-default}

@parblock
```cpp
template <typename T1, typename T2>
void dot(const T1& a, const T2& b,  Traits<T1>::scalar_type& result)

template <typename T1, typename T2>
Traits<T1>::scalar_type dot(const T1& a, const T2& b)
```
- Computes the dot product of `a` and `b` storing the result in `result` (first variant) or returning it directly (second variant).
@endparblock

#### elementwise_multiply()

@m_class{m-block m-default}

@parblock
```cpp
template <typename T, typename T1, typename T2>
void elementwise_multiply(
  Traits<T>::scalar_type alpha,
  const T& x,
  const T1& z,
  Traits<T>::scalar_type beta,
  T2& y)
 ```
- Multiplies `x` and `z` vectors element-wise computing: @f$y_i := \alpha x_i z_i + \beta y_i@f$
@endparblock

#### update()

@m_class{m-block m-default}

@parblock
```cpp
template<class T, class T1, class scalar_t>
void update(T& v, scalar_t alpha,
  const T1& u1, scalar_t b1)

template<class T, class T1, class T2, class scalar_t>
void update(T& v, scalar_t alpha,
  const T1& u1, scalar_t b1,
  const T2& u2, scalar_t b2)

template<class T, class T1, class T2, class T3, class scalar_t>
void update(T& v, scalar_t alpha,
  const T1& u1, scalar_t b1,
  const T2& u2, scalar_t b2,
  const T3& u3, scalar_t b3)

template<class T, class T1, class T2, class T3,  class T4, class scalar_t>
void update(T& v, scalar_t alpha,
  const T1& u1, scalar_t b1,
  const T2& u2, scalar_t b2,
  const T3& u3, scalar_t b3,
  const T4& u4, scalar_t b4)
```
- Updates `v` with `u1`-`u4` (all vectors or matrices) computing @f$v_\lambda := \alpha v_\lambda + \sum_{j=1}^k{\beta_j u_{j\lambda}}@f$<br> where @f$\lambda@f$ is rank one or rank two index and @f$k=1,2,3,4@f$ depends on the selected overload.
@endparblock

#### product()

@m_class{m-block m-default}

@parblock
```cpp
template <typename Mode, typename scalar_type,
  class A_type, class x_type, class y_type>
void product(
  Mode /* mode */,
  const scalar_type alpha,
  const A_type& A,
  const x_type& x,
  const scalar_type beta,
  y_type& y)
 ```
- Multiplies matrix `A` and vector `x` computing @f$y := \alpha A x + \beta y@f$ for `Mode = pressio::nontranspose` and @f$y := \alpha A^T x + \beta y@f$ for `Mode = pressio::transpose`.
@endparblock

@m_class{m-block m-default}

@parblock
```cpp
template <typename ModeA, typename ModeB, typename scalar_type,
  class A_type, class B_type, class C_type>
void product(
  ModeA /* modeA */,
  ModeB /* modeB */,
  const scalar_type alpha,
  const A_type& A,
  const B_type& B,
  const scalar_type beta,
  C_type& C)
```
- Multiplies matrices `A` and `B` computing @f$C := \alpha A B + \beta C@f$ if both modes are `pressio::nontranspose`. Use `pressio::transpose` to transpose input matrices of choice.
@endparblock

@m_class{m-block m-default}

@parblock
```cpp
template <class scalar_type, class A_type, class C_type>
void product(
  ::pressio::transpose /* modeA */,
  ::pressio::nontranspose /* modeB */,
  const scalar_type alpha,
  const A_type& A,
  const scalar_type beta,
  C_type& C)
```
- Computes @f$C := \alpha A^T A + \beta C@f$
@endparblock

@m_class{m-block m-default}

@parblock
```cpp
template <class C_type, class scalar_type, class A_type>
C_type product(
  ::pressio::transpose /* modeA */,
  ::pressio::nontranspose /* modeB */,
  const scalar_type alpha,
  const A_type& A)
```
- Computes and returns @f$C = \alpha A^T A@f$
@endparblock


## Utilities

#### extent()

@m_class{m-block m-default}

@parblock
```cpp
template<class T, class IndexType>
Traits<T>::size_type extent(const T& objectIn, const IndexType i)
```
- Returns the number of elements that `objectIn` has in `i`-th dimension (zero based).
@endparblock

#### matching_extents::compare()

@m_class{m-block m-default}

@parblock
```cpp
template<typename T1, typename T2>
bool matching_extents<T1, T2>::compare(const T1& a, const T2& b)
```
- Returns `true` if objects `a` and `b` have same number of elements in each dimension.
@endparblock

#### clone()

@m_class{m-block m-default}

@parblock
```cpp
template <class T>
T clone(const T& src)
```
- Creates a copy of `src` object.
@endparblock

#### deep_copy()

@m_class{m-block m-default}

@parblock
```cpp
template<typename T1, typename T2>
void deep_copy(T1& dest, const T2& src)
```
- Copies `src` to `dest` assuring a deep copy.
@endparblock

#### fill() and set_zero()

@m_class{m-block m-default}

@parblock
```cpp
template <typename T>
void fill(T& o, Traits<T>::scalar_type value)

template <typename T>
void set_zero(T& o)
```
- `fill()` populates `o` (matrix or vector) with given scalar `value`.
- `set_zero()` populates `o` (matrix or vector) with zeros.
@endparblock

@m_class{m-block m-default}

#### resize()

@m_class{m-block m-default}

@parblock
```cpp
template <typename T>
void resize(T& o, const Traits<T>::size_type newSize)
```
- Resizes vector `o` into `newSize` number of elements.
@endparblock

@m_class{m-block m-default}

@parblock
```cpp
template <typename T>
void resize(
  T& o,
  const Traits<T>::size_type newRows,
  const Traits<T>::size_type newCols)
```
- Resizes matrix `o` into `newCols` x `newRows` number of elements.
@endparblock
