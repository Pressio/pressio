
# Summary

Description of every test problem in this directory.
If you add one, make sure you update this file.

References for useful problems are:

- http://www2.imm.dtu.dk/documents/ftp/tr99/tr05_99.pdf

- https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml

- http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf

For each problem, we provide (not necessarily in this order):
description, link to header, data types used, which API meets,
why it is included, where it was found

========================================================
# problem 1
========================================================
- determined system defined as:
  ```cpp
  r(0) = x^3 + y - 1.
  r(1) = -x + y^3 + 1.
  ```
  with exact jacobian provided

- [header](./problem1.hpp)

- Eigen

- RealSystemWithFusedResidualAndJacobian

- very simple problem

- made up or found somewhere random

========================================================
# problem 2
========================================================
- determined system, intersection of a circle and hyperbola,
  but with an inexact jacobian

- [header](./problem2.hpp)

- Eigen

- RealSystemWithResidualAndJacobian

- inexact jacobian

- pg.128 here http://www.math.iit.edu/~fass/477577_Chapter_17.pdf

========================================================
# problem 3
========================================================
- data fit problem with 2 unknowns and 8 residuals (overdetermined)

- [header](./problem3.hpp)

- Eigen

- RealSystemWithResidualAndJacobian

- interesting data fitting problem

- not sure where it is from

========================================================
# problem 4
========================================================
- Rosenbrock with 3 vars (overdetermined)

- [header](./problem4.hpp)

- Eigen

- RealSystemWithResidualAndJacobian

- interesting because Rosenbrock function is widely used

- can find on wikipedia or somewhere else

========================================================
# problem 5
========================================================
- Rosenbrock with 4 vars (overdetermined)

- [header](./problem5.hpp)

- Eigen

- APIs:
  - Problem5a: RealSystemWithResidualAndJacobian

- interesting because Rosenbrock function is widely used

- can find on wikipedia or somewhere else

========================================================
# problem 6
========================================================
- datafit with 5 vars and 33 equations (overdetermined)
- [header](./problem6.hpp)
- Eigen
- RealSystemWithResidualAndJacobian
- minpack problem
- links:
   - sec 3.4 http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
   - data from http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dedffj.f

========================================================
# problem 7
========================================================
- datafit with 5 vars and 33 equations (overdetermined)
- [header](./problem6.hpp)
- Tpetra
- RealSystemWithResidualAndJacobian
- minpack problem
- links:
   - sec 3.4 http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
   - data from http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dedffj.f


========================================================
# problem 8
========================================================
- nist problem (overdetermined)

- [header](./problem8.hpp)

- Eigen

- RealSystemWithResidualAndJacobian

- hard problem

- mgh10, from: https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
   see also: https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/MGH10.dat

========================================================
# problem 9
========================================================
- datafit with 11 vars and 65 equations (overdetermined)

- [header](./problem7.hpp)

- Eigen

- RealSystemWithResidualAndJacobian

- minpack problem

- links:
   - sec 3.5 http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
     (note that there is a typo: n=9 is wrong, should be n=11
   - data from http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dgdffj.f
