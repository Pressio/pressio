
# mpl


@m_class{m-note m-default}

@parblock
Defined in header: `<pressio/mpl.hpp>`

Public namespace: `pressio::mpl`
@endparblock

## Overview

@m_class{m-note m-info}

@parblock
Provides metaprogramming functionalities that are useful for a variety of purposes
and are one of the fundamental building blocks for the other parts of the pressio library.
@endparblock

If you are familiar with the `<type_traits>` header from
the standard library, the `pressio/mpl` will look familiar too.
Some parts of `pressio/mpl` have been adapted from the [tinympl project](http://sbabbi.github.io/tinympl).
The tinympl project appears to be no longer maintained.


## Content

The following is a *partial* list only intended to provide a general idea of the supported features.

To find out all supported cases, browse the [source](https://github.com/Pressio/pressio/tree/main/include/pressio/mpl),
and for example usage, see [this](https://github.com/Pressio/pressio/blob/main/tests/functional_small/mpl/all.cc).


@m_class{m-block m-default}

@parblock
```cpp
template<class T> struct not_void;
```

- Checks if a type `T` is NOT a void type.
  Provides the static member constant `value` that is equal to true, if `T` is NOT of
  the type `void`, `const void`, `volatile void`, or `const volatile void`.
  Otherwise, value is true.

- Example:<br/>
  ```cpp
  namespace pmpl = pressio::mpl;
  static_assert(pmpl::not_void<double>::value, "" );
  ```
@endparblock


@m_class{m-block m-default}

@parblock
```cpp
template< template<class ... T> class F, class ... Args> struct all_of;
```

- Determines whether every element in the sequence satisfies the given predicate.
  The predicate `F` must be such that `F<T>::value` must be convertible to `bool`.
  Provides the static member constant `value` that is equal to true iff
  all the elements in the sequence satisfy the predicate `F`.
  Otherwise, value is false.

- Example:<br/>
  ```cpp
  namespace pmpl = pressio::mpl;
  static_assert(pmpl::all_of<std::is_floating_point, double, float>::value, "" );
  ```
@endparblock


@m_class{m-block m-default}

@parblock
```cpp
template< template<class ... T> class F, class ... Args> struct any_of;
```

- Determines whether any element in the sequence satisfies the given predicate.
  The predicate `F` must be such that `F<T>::value` must be convertible to `bool`.
  Provides the static member constant `value` that is equal to true iff
  at least one element in the sequence satisfies the predicate `F`.
  Otherwise, value is equal to false.
@endparblock


@m_class{m-block m-default}

@parblock
```cpp
template< template<class ... T> class F, class ... Args> struct none_of;
```

- Determines whether none of the elements in the sequence satisfy the given predicate.
  The predicate `F` must be such that `F<T>::value` must be convertible to `bool`.
  Provides the static member constant `value` that is equal to true iff
  none of the elements in the sequence satisfy the predicate `F`.
  Otherwise, value is equal to false.
@endparblock


@m_class{m-block m-default}

@parblock
```cpp
template<class T, class IndexType> struct is_subscriptable_as;
```

- Provides the static member constant `value` that is equal to true if
  `T` has subscript operator `[]`, it can be indexed by an instance of `IndexType`,
  and the return type is not void.
  Otherwise, value is equal to false.
@endparblock
