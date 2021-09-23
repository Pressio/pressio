
# utils

@m_class{m-note m-default}

@parblock
Defined in header: `<pressio/utils.hpp>`

Public namespaces: `pressio`, `pressio::utils`
@endparblock


## Logger

One of the main functionalities inside `utils` is the logger.
To implement the pressio logging functionalities, we have used [spdlog](https://github.com/gabime/spdlog),
modified a subset of it, for example to work seamlessly with MPI.

@m_class{m-note m-warning}

@parblock
By default, for performance reasons, the logger is disabled, so pressio will not output anything.
@endparblock

<!-- If you just `#include<pressio/utils.hpp>` and expect logging, you will be disappointed!<br/> -->
<!-- @m_span{m-text m-warning}By default, for performance reasons, the logger is disabled.@m_endspan -->

@m_class{m-note m-info}

@parblock
To enable logging, you need to do two things:
1. insert a `define` statement to set the *minimum* level *before* including the utils header
2. initialize and finalize the logger singleton
@endparblock

The following snippet provides the main idea:

```cpp
// this sets the default min level
#define PRESSIO_LOG_ACTIVE_MIN_LEVEL 2

// this is needed since it contains the logging code
#include <pressio/utils.hpp>

int main()
{
  // namespace alias for readability
  namespace plog = pressio::log;

  // initialize logger (see below for details)

  // your code, use logger

  // finalize logger (see below for details)
}
```

### Initializing and finalizing

To initialize, you have these choices:

```cpp
// only log messages to terminal
plog::initialize(pressio::logto::terminal);

// log messages to terminal and a file named, e.g., "my_log.txt"
plog::initialize(pressio::logto::fileAndTerminal, "my_log.txt");

// only log messages to a file named, e.g., "my_log.txt"
plog::initialize(pressio::logto::file, "my_log.txt");
```

And to finalize:

```cpp
plog::finalize();
```

### If running with MPI

If you are running with MPI, the logger prints to the terminal *only from rank==0*.
However, it automatically creates a per-rank log file if you choose the file output.
For example, the following code:

```cpp
int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int my_rank = {};
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  namespace plog = pressio::log;
  plog::initialize(pressio::logto::file, "log_file.txt");
  PRESSIOLOG_INFO("print from rank {:2}", my_rank);
  plog::finalize();

  MPI_Finalize();
}
```

If we were to run this with N ranks, we would obtain two
files `log_file.txt_0`, and `log_file.txt_1`.
Currently, the logger works only for the world communicator.
We will later extend the API to accept a communicator object.




### Levels

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
but then at runtime you can reset for a *higher* level (see below).
Note that you cannot reset the level to something that is *lower* than the
one you set via the `define` statement.

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
  // initialize logger

  double a = 1.1;
  PRESSIOLOG_TRACE("my value is {:.6f}", a);
  PRESSIOLOG_DEBUG("my value is {:.6f}", a);
  PRESSIOLOG_INFO("my value is {:.6f}", a);
  PRESSIOLOG_WARN("my value is {:.6f}", a);
  PRESSIOLOG_ERROR("my value is {:.6f}", a);
  PRESSIOLOG_CRITICAL("my value is {:.6f}", 55.6);

  // finalize logger
}
```

where we note that you can use the [{fmt} library](https://github.com/fmtlib/fmt)
to properly format the print statements.

@m_class{m-block m-warning}

@par Keep in mind:
The log statements issued for a specific level will be printed
*only if* `PRESSIO_LOG_ACTIVE_MIN_LEVEL` is smaller or equal than that level.
If the logger is disabled, the macros are expanded to a no-op.
So it does not cost you anything to place log statements in your code,
because in production mode you can just compile to no-op.
