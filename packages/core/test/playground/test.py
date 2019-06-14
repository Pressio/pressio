
import main
import numpy as np

class MyFom:
  def __init__(self, params):
      self.params = params

  def residual(self, A, t):
      R = np.zeros(5)
      print( "fromPy", len(A) )
      print(A[0][0])
      A[0][0] = 1.1
      return R

appObj = MyFom("dummy")
main.test(appObj)
#print(a)

main.foo(np.ones(3))
# s = np.zeros(5)
# s[:] = 1
# print(s)

# appObj = MyFom("dummy")
# r1 = appObj.residual(s, 0)
# print(r1)
