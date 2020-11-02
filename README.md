
# Overview

Pressio is a collection of repositories providing reduced-order models (ROMs) capabilties.
It is a open-source project aimed at enabling leading-edge projection-based
reduced order models (pROMs) for dynamical systems in science and engineering.

# Userguide and documentation
This repository is specific for the C++ library.
You can find the full user-guide [here](https://pressio.github.io/pressio/html/index.html).

The full Pressio project site can be found [here](https://pressio.github.io).

# License and Citation
The full license is available [here](https://pressio.github.io/various/license/).

We are working on publishing this: you can find our arXiv preprint at: https://arxiv.org/abs/2003.07798


<!-- * `pressio-tutorials`: C++ tutorials explaining how to use `pressio` and its functionalities; -->

<!-- * `pressio-builder`: an auxiliary repo with bash helper scripts for configuring/building/installing `pressio`, and `pressio-tutorials`. -->

<!-- ## Questions -->
<!-- For questions, find us on Slack: https://pressioteam.slack.com or open an issue. -->

<!-- ## License and Citation -->
<!-- Pressio is released with the following [LICENSE](./LICENSE). -->

<!-- Please see the following axXiv paper: https://arxiv.org/abs/2003.07798 -->

<!-- ## Structure -->
<!-- For a description of `pressio` code structure, see [here](https://github.com/Pressio/pressio/wiki/Structure-of-pressio). -->

<!-- ## Building and Installing -->

<!-- ### If you only want to use `pressio` from your code -->
<!-- In this case, since `pressio` is header-only, there is **no building process needed**. -->
<!-- You clone the `pressio` repo, and within your code you include the `pressio/packages` to find the `pressio` headers. -->
<!-- However, since `pressio` uses preprocessor directives to selectively enable/disable code for target TPLs, when you build your code you need to have these preprocessor directives defined. -->
<!-- For example, if your code uses Trilinos, to enabled the Trilinos-related code in `pressio` you need to have `PRESSIO_ENABLE_TPL_TRILINOS` defined *before* you include -->
<!-- the `pressio` headers. The list of CMake options to enable can be found [here](./list_of_cmake_optional_vars_to_enable.md). -->

<!-- ### If you want to build the unit and regression tests in `pressio` -->
<!-- Sample cmake configure lines can be found [here](https://github.com/Pressio/pressio/wiki/Sample-CMake-configure-lines-for-pressio). -->

<!-- Follow [this](https://github.com/Pressio/pressio/wiki/Serial-build-of-Pressio-with-tests-enabled) for a basic *serial* build that uses only GTest and Eigen and it is done with `pressio-builder` (which automatically builds) Gtest, Eigen for you. -->

<!-- ## Sample codes -->
<!-- While we improve the tutorials, please look at the subdirectory `pressio/tests` -->

<!-- ## Disclaimer -->

<!-- * Pressio is work-in-progress. At the time of this writing, it is a fairly young project and things are obviously evolving. Several package would benefit from substantial work on testing and documentation, and this is ongoing. However, `pressio` is functional and has been already tested/deployed on large-scale applications. -->
