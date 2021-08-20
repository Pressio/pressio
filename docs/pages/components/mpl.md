
# mpl

Defined in header: `<pressio/mpl.hpp>`

Public namespace: `pressio::mpl`

## Overview

Provides metaprogramming functionalities that are useful for a variety of purposes
and are one of the fundamental building blocks for the other parts of the pressio library.
If you are familiar with the `<type_traits>` header from
the standard library, the `pressio/mpl` will look familiar too.

Some parts of `pressio/mpl` have been adapted/extended/forked
from the tinympl project, forked from http://sbabbi.github.io/tinympl.
The tinympl project appears to be no longer maintained.

## Content

We are not going to list all the metafunctions supported, but we only
give a glimpse of some of them.

### pressio::mpl::not_void

```cpp
template<class T> struct not_void;
```

Checks whether a type `T` is NOT a void type.
Provides the static member constant `value` that is equal to true, if `T` is NOT of
the type `void`, `const void`, `volatile void`, or `const volatile void`.
Otherwise, value is equal to true.


### pressio::mpl::is_subscriptable_as

```cpp
template<class T, class IndexType> struct is_subscriptable_as;
```

Provides the static member constant `value` that is equal to true if
`T` has subscript operator `[]`, it can be indexed by an instance of `IndexType`,
and the return type is not void.
Otherwise, value is equal to false.


### to do: add more
