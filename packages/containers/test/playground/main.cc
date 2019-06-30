
#include <iostream>
#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <functional>
#include <pybind11/eigen.h>

namespace py = pybind11;

template <typename T>
void foo(T a){
  std::cout << "foo_py " << a(0) << std::endl;
}


//void run_test(const std::function<int(int)>& f)
void test(py::object & f)
{
  Eigen::MatrixXd A(2,2);
  A.setConstant(2.2);

  const py::object & f_r = f.attr("residual");
  f_r(A, 0);
  std::cout << "GG " << A(0,0) << std::endl;
  //py::print( f_r(0) );
  //std::cout << "result: " << r << std::endl;
}

PYBIND11_MODULE(main, m) {
  m.def("test", &test, "adds two numbers");

  m.def("foo", &foo<Eigen::MatrixXd>, "dkdk");
}

// // int add(int i, int j) {
// //     return i + j;
// // }

// PYBIND11_MODULE(main, m) {
//     m.doc() = "pybind11 example plugin";
//     m.def("add", &add, "A function which adds two numbers");
// }
