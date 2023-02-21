
# Summary

Description of every test problem in this directory.
If you add one, make sure you update this file.


Some references for useful problems are:

- http://www2.imm.dtu.dk/documents/ftp/tr99/tr05_99.pdf

- https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml

- https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml


Legend:
1. link to header and data types
2. a description of the problem
3. which API meets
4. why it is included
5. where it was found

========================================================
# problem 1

1. [header](./problem1.hpp), Eigen

2. description:
  ```
  r(0) = x^3 + y - 1.
  r(1) = -x + y^3 + 1.
  ```
  and exact jacobian is provided

3. `DeterminedSystemWithRealRhsAndJacobian` and
   `DeterminedSystemWithFusedRealRhsAndJacobian`

4. it is simple

5. don't remember

========================================================
# problem 2

1. [header](./problem2.hpp), Eigen

2. intersection of a circle and hyperbola, but with an inexact jacobian

3. `DeterminedSystemWithRealRhsAndJacobian`

4. inexact jacobian

5. pg.128 here http://www.math.iit.edu/~fass/477577_Chapter_17.pdf

========================================================
# problem 3

1. [header](./problem3.hpp), Eigen

2. data fit problem with 2 unknowns and 8 residuals

3. `OverDeterminedSystemWithRealRhsAndJacobian`

4. interesting data fitting problem

5. not sure where it is from

========================================================
# problem 4

1. [header](./problem4.hpp), Eigen

2. Rosenbrock with 3 vars

3. `OverDeterminedSystemWithRealRhsAndJacobian`

4. interesting because Rosenbrock function is widely used

5. can find on wikipedia or somewhere else

========================================================
# problem 5

1. [header](./problem5.hpp), Eigen

2. Rosenbrock with 4 vars

3. we have:
   - `Problem5a: `OverDeterminedSystemWithRealRhsAndJacobian`

4. interesting because Rosenbrock function is widely used

5. can find on wikipedia or somewhere else

========================================================
# problem 6

1. [header](./problem6.hpp), Epetra

2. datafit with 5 vars and 33 equations

3. `Problem6: `OverDeterminedSystemWithRealRhsAndJacobian`

4. interesting because it is minpack problem

5. links:
   - sec 3.4 http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
   - data from http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dedffj.f

========================================================
# problem 7

1. [header](./problem7.hpp), Epetra

2. datafit with 11 vars and 65 equations

3. `Problem7: `OverDeterminedSystemWithRealRhsAndJacobian`

4. interesting because it is minpack problem

5. links:
   - sec 3.5 http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
   - data from http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dgdffj.f

========================================================
# problem 8

1. [header](./problem8.hpp), Eigen

2. nist problem

3. `OverDeterminedSystemWithRealRhsAndJacobian`

4. interesting because it is a hard problem

5. mgh10, from: https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
   see also: https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/MGH10.dat
