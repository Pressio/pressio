## Docker images
Docker images which could be found in docker_scripts directory are used in GitHub Actions for testing purposes

### Naming convention
Let's decode an image name:

`ubuntu-20.04-gnu_9-eigen_3.3.7-gtest-trilinos_9fec35276d846a667bc668ff4cbdfd8be0dfea08`

Decoded:
- `ubuntu-20.04` - base image and OS
- `gnu_9` - compilers used in the image are from GNU Compiler Collection version 9 in this case
- `eigen_3.3.7` - Eigen is C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms in this case version 3.3.7
- `gtest` - GoogleTest, newest version available from master branch on dockerimage creation is used
- `trilinos_9fec35276d846a667bc668ff4cbdfd8be0dfea08` - Trilinos version from `9fec35276d846a667bc668ff4cbdfd8be0dfea08` hash is build in this case that corresponds to Trilinos release 13.0.0

Intel compilers image:
`intel_oneapi-eigen_3.3.7-gtest`

Decoded:
- `intel_oneapi` - image is based on `intel/oneapi-hpckit` which is based on Ubuntu 18.04
- rest is exactly the same as in example decoded above

### Usage
Images described above are used for job named `CI`.

Images:
- `ubuntu-20.04-gnu_10-eigen_3.3.7-gtest-trilinos_rel_13.0.0`
- `ubuntu-20.04-gnu_10-eigen_3.3.7-gtest-trilinos_rel_13.0.1`

are used for job named `CI-trilinos`, which shall be soon defined differently. Naming convention is the same, but at the end there is a Trilinos release.