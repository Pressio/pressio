.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst


``mpl``
=======

.. admonition:: Info
   :class: important

   Header: ``<pressio/mpl.hpp>``

   Public namespace: ``pressio::mpl``


.. admonition:: Description

    :medium:`Provides metaprogramming functionalities that are always useful
    for generic programming and are one of the fundamental
    building blocks of the pressio library.`


If you are familiar with the ``<type_traits>`` header from
the standard library, the ``pressio/mpl`` will look familiar too.
Some parts of ``pressio/mpl`` have been adapted from the `tinympl project <http://sbabbi.github.io/tinympl>`_.
The tinympl project appears to be no longer maintained.

Content
-------

The following is a *partial* list only intended to provide a general idea of the supported features.

To find out all supported cases, browse the `source <https://github.com/Pressio/pressio/tree/main/include/pressio/mpl>`__.


``not_void``
~~~~~~~~~~~~

.. code-block:: cpp

   template<class T> struct not_void;

*
  Provides the static member constant ``value`` that is equal to true, if ``T`` is NOT of
  the type ``void``\ , ``const void``\ , ``volatile void``\ , or ``const volatile void``.
  Otherwise, value is true.

*
  Example:\ :raw-html-m2r:`<br/>`

  .. code-block:: cpp

     namespace pmpl = pressio::mpl;
     static_assert(pmpl::not_void<double>::value, "" );


``all_of``
~~~~~~~~~~

.. code-block:: cpp

   template< template<class ... T> class F, class ... Args> struct all_of;

*
  Determines whether every element in the sequence satisfies the given predicate.
  The predicate ``F`` must be such that ``F<T>::value`` must be convertible to ``bool``.
  Provides the static member constant ``value`` that is equal to true iff
  all the elements in the sequence satisfy the predicate ``F``.
  Otherwise, value is false.

*
  Example:\ :raw-html-m2r:`<br/>`

  .. code-block:: cpp

     namespace pmpl = pressio::mpl;
     static_assert(pmpl::all_of<std::is_floating_point, double, float>::value, "" );


``any_of``
~~~~~~~~~~

.. code-block:: cpp

   template< template<class ... T> class F, class ... Args> struct any_of;

* Determines whether any element in the sequence satisfies the given predicate.
  The predicate ``F`` must be such that ``F<T>::value`` must be convertible to ``bool``.
  Provides the static member constant ``value`` that is equal to true iff
  at least one element in the sequence satisfies the predicate ``F``.
  Otherwise, value is equal to false.


``none_of``
~~~~~~~~~~~

.. code-block:: cpp

   template< template<class ... T> class F, class ... Args> struct none_of;

* Determines whether none of the elements in the sequence satisfy the given predicate.
  The predicate ``F`` must be such that ``F<T>::value`` must be convertible to ``bool``.
  Provides the static member constant ``value`` that is equal to true iff
  none of the elements in the sequence satisfy the predicate ``F``.
  Otherwise, value is equal to false.


``is_subscriptable_as``
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   template<class T, class IndexType> struct is_subscriptable_as;

* Provides the static member constant ``value`` that is equal to true if
  ``T`` has subscript operator ``[]``\ , it can be indexed by an instance of ``IndexType``\ ,
  and the return type is not void.
  Otherwise, value is equal to false.
