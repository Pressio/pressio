
# utils

Defined in header: `<pressio/utils.hpp>`

Public namespaces: `pressio` and `pressio::utils`.


## Logger

One of the main functionalities inside `utils` is the logger.
To implement the pressio logging functionalities, we have used [spdlog](https://github.com/gabime/spdlog),
and rewritten a subset of it, for example to work seamlessly with MPI.

@m_class{m-block m-warning}

@par By default all logging is disabled.
If you just `#include<pressio/utils.hpp>` and expect logging, you will be disappointed!
By default, all logging is disabled for performance reasons.
To enable and use it, you need to do two things:
1. place a `define` statement to set the desired level: note that
this define statement *must* be placed *before* including the utils header.

2. you need to start and finalize the logger singleton.

As follows:
```cpp
// this sets the default min level
#define PRESSIO_LOG_ACTIVE_MIN_LEVEL 2

// this is needed to include the logger code
#include <pressio/utils.hpp>

int main()
{
  namespace plog   = pressio::log;
  plog::initialize(pressio::logto::terminal);

  // your code

  plog::finalize();
}
```
In this case, we set the level to 2, which means "info"-level logging.
The supported levels are:

```cpp
#define PRESSIO_LOG_LEVEL_TRACE		0
#define PRESSIO_LOG_LEVEL_DEBUG		1
#define PRESSIO_LOG_LEVEL_INFO		2
#define PRESSIO_LOG_LEVEL_WARN		3
#define PRESSIO_LOG_LEVEL_ERROR		4
#define PRESSIO_LOG_LEVEL_CRITICAL	5
#define PRESSIO_LOG_LEVEL_OFF		6
```

### Resetting the level

If you want, you can use the define statement to set the min level,
but then are runtime reset it as follows:
```cpp
  // your code
  // ...
  plog::setVerbosity({plog::level::info});
  // ...
```


### The loggin macros

To actually issue log statements, you use the macros as in the following example:
```cpp
int main()
{
  double a = 1.1;
  PRESSIOLOG_TRACE("my value is", a, ");
  PRESSIOLOG_DEBUG("my value is", a, " but I am debug");
  PRESSIOLOG_INFO("my value is", a, " but I am info");
  PRESSIOLOG_WARN("my value is", a, " but I am warn");
  PRESSIOLOG_ERROR("my value is", a, " but I am an error level");
  PRESSIOLOG_CRITICAL("my value is", 55.6, " but don't forget );
}
```

The log statements issued for a specific level will be printed *only if* `PRESSIO_LOG_ACTIVE_MIN_LEVEL`
is greater or equal than that level.
