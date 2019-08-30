
#Pressio Coding Style

##Indentation
2 tabs


##Formatting
- the width of the code should not be exceeding about 75 characters.
  Favor vertical development rather than horizontal one.

- Prefer wrapping text around than huge long lines, or it will be hard 
  for other developers to read code because everybody uses different editors.

- use some spacing where appropriate for clarity. This helps readability and to keep things close together that actually belong together. When you have spaces between lines, please use 1 row spacing. Usually, >=2 lines of spacing is just too much. 

- preferably put comments above rather than on the right-hand-side of statements. 
E.g., prefer: 
```C++
// set tolerance
const double tol = 1e-6;
```
rather than: 
```C++
const double tol = 1e-6; // set tolerance
```


##Filenames
- A file name should conform to:   packageName_name_describing_purpose_of_file
- Prepend every file with the name of the package that file belongs to, followed by the specific name of the class or whatever that file covers.

- file names should be expressive. From the filename I should be able to guess what it contains.

- If the name of the class is made of multiple words, then the file name should be all lowercase with underscore separating the components words making up the class but the class name itself should be ONE single string of words where each subword starts with uppercare letter. 
For instance, if the target class name is: 
```
VectorEpetraDistributed
```
and the file belongs to the package "core", then the filename should be:
```
core_vector_epetra_distributed.hpp
```


##Classes
- class names should beging with uppercase letter
- no spaces in the names, user camelCase convention
- methods should beging with lowercase letter


##The #define Guard
All header files should have #define guards. The format of the symbol name should be:
```
<PACKAGENAME>_<FILE>_HPP_
```
For example, the file 
```
packages/core/src/vector/core_vector_traits.hpp
```
should have the following guard:
```
#ifndef CORE_VECTOR_TRAITS_HPP_
#define CORE_VECTOR_TRAITS_HPP_

...

#endif  //CORE_VECTOR_TRAITS_HPP_
```

##Namespaces
- Namespaces should have unique names and should be lowercase.
- Do not use using-directives (e.g. using namespace blabla), unless it is confined within a well-defined scope that is not going to pollute the overall project.
- after definition of namespace, no indentation. For example:
```cpp
namespace test{
  
struct Foo{}
...
template<typename scalar>
class Vector{
	...
};

void fooFunction(){
	...
}
  
}
```

##File organization
- The src directory of a package should not have a bunch of files all dumped there. 
Things should be separated logically based on functionality. 
- Best practice: create subfolders to contain classes that belong together. 
- The CMakeList inside the src folder should preferably be written using globs to make the process of adding headers much easier and avoiding forgetting files. 




