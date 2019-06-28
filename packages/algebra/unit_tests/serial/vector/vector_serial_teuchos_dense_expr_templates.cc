
#include <gtest/gtest.h>
#include "ALGEBRA_VECTOR"
#include "ALGEBRA_OPS"

using natvec_t = Teuchos::SerialDenseVector<int, double>;
using myvec_t = ::rompp::algebra::Vector<natvec_t>;

TEST(algebra_vector_teuchos_serial_dense_class,
     additionOperator){

  myvec_t v1(4);
  v1[0] = 3.; v1[1] = 2.;
  v1[2] = 4.; v1[3] = 5.;
  myvec_t v2(4);
  v2[0] = 1.; v2[1] = 1.;
  v2[2] = 1.; v2[3] = 1.;

  myvec_t v3(4);
  v3 = (v1+v2);
  EXPECT_DOUBLE_EQ( v3[0], 4.0);
  EXPECT_DOUBLE_EQ( v3[1], 3.0);
  EXPECT_DOUBLE_EQ( v3[2], 5.0);
  EXPECT_DOUBLE_EQ( v3[3], 6.0);

  myvec_t v4( v1+v2 );
  EXPECT_DOUBLE_EQ( v4[0], 4.0);
  EXPECT_DOUBLE_EQ( v4[1], 3.0);
  EXPECT_DOUBLE_EQ( v4[2], 5.0);
  EXPECT_DOUBLE_EQ( v4[3], 6.0);
}


TEST(algebra_vector_teuchos_serial_dense_class,
     subtractOperator){

  myvec_t v1(4);
  v1[0] = 3.; v1[1] = 2.;
  v1[2] = 4.; v1[3] = 5.;
  myvec_t v2(4);
  v2[0] = 1.; v2[1] = 1.;
  v2[2] = 1.; v2[3] = 1.;

  myvec_t v3(4);
  v3 = (v1-v2);
  EXPECT_DOUBLE_EQ( v3[0], 2.0);
  EXPECT_DOUBLE_EQ( v3[1], 1.0);
  EXPECT_DOUBLE_EQ( v3[2], 3.0);
  EXPECT_DOUBLE_EQ( v3[3], 4.0);

  myvec_t v4( v1-v2 );
  EXPECT_DOUBLE_EQ( v4[0], 2.0);
  EXPECT_DOUBLE_EQ( v4[1], 1.0);
  EXPECT_DOUBLE_EQ( v4[2], 3.0);
  EXPECT_DOUBLE_EQ( v4[3], 4.0);
}


TEST(algebra_vector_teuchos_serial_dense_class,
     timesOperatorWithVector){

  myvec_t v1(4);
  v1[0] = 3.; v1[1] = 2.;
  v1[2] = 4.; v1[3] = 5.;
  myvec_t v2(4);
  v2[0] = 1.; v2[1] = 1.;
  v2[2] = 1.; v2[3] = 1.;

  myvec_t v3(4);
  v3 = v1*v2;
  EXPECT_DOUBLE_EQ( v3[0], 3.0);
  EXPECT_DOUBLE_EQ( v3[1], 2.0);
  EXPECT_DOUBLE_EQ( v3[2], 4.0);
  EXPECT_DOUBLE_EQ( v3[3], 5.0);

  myvec_t v4( v1*v2 );
  EXPECT_DOUBLE_EQ( v4[0], 3.0);
  EXPECT_DOUBLE_EQ( v4[1], 2.0);
  EXPECT_DOUBLE_EQ( v4[2], 4.0);
  EXPECT_DOUBLE_EQ( v4[3], 5.0);
}

TEST(algebra_vector_teuchos_serial_dense_class,
     timesOperator3Vectors){

  myvec_t v1(4);
  v1[0] = 3.; v1[1] = 2.;
  v1[2] = 4.; v1[3] = 5.;
  myvec_t v2(4);
  v2[0] = 1.; v2[1] = 1.;
  v2[2] = 1.; v2[3] = 1.;
  myvec_t v3(4);
  v3[0] = 2.; v3[1] = 1.;
  v3[2] = 1.; v3[3] = 3.;

  myvec_t v4(4);
  v4 = v1*v2*v3;
  EXPECT_DOUBLE_EQ( v4[0], 6.0);
  EXPECT_DOUBLE_EQ( v4[1], 2.0);
  EXPECT_DOUBLE_EQ( v4[2], 4.0);
  EXPECT_DOUBLE_EQ( v4[3], 15.0);

  myvec_t v5( v1*v2*v3 );
  EXPECT_DOUBLE_EQ( v5[0], 6.0);
  EXPECT_DOUBLE_EQ( v5[1], 2.0);
  EXPECT_DOUBLE_EQ( v5[2], 4.0);
  EXPECT_DOUBLE_EQ( v5[3], 15.0);
}


TEST(algebra_vector_teuchos_serial_dense_class,
     timesOperatorWithScalar){

  myvec_t v1(4);
  v1[0] = 3.; v1[1] = 2.;
  v1[2] = 4.; v1[3] = 5.;

  constexpr double b = 2.;

  myvec_t v3(4);
  v3 = v1*b;
  EXPECT_DOUBLE_EQ( v3[0], 6.0);
  EXPECT_DOUBLE_EQ( v3[1], 4.0);
  EXPECT_DOUBLE_EQ( v3[2], 8.0);
  EXPECT_DOUBLE_EQ( v3[3], 10.0);

  myvec_t v4(4);
  v4 = (b*v1);
  EXPECT_DOUBLE_EQ( v4[0], 6.0);
  EXPECT_DOUBLE_EQ( v4[1], 4.0);
  EXPECT_DOUBLE_EQ( v4[2], 8.0);
  EXPECT_DOUBLE_EQ( v4[3], 10.0);

  myvec_t v5( v1*b );
  EXPECT_DOUBLE_EQ( v5[0], 6.0);
  EXPECT_DOUBLE_EQ( v5[1], 4.0);
  EXPECT_DOUBLE_EQ( v5[2], 8.0);
  EXPECT_DOUBLE_EQ( v5[3], 10.0);

  myvec_t v6( b*v1);
  EXPECT_DOUBLE_EQ( v6[0], 6.0);
  EXPECT_DOUBLE_EQ( v6[1], 4.0);
  EXPECT_DOUBLE_EQ( v6[2], 8.0);
  EXPECT_DOUBLE_EQ( v6[3], 10.0);
}


TEST(algebra_vector_teuchos_serial_dense_class,
     compoundAssignPlus){

  myvec_t a(4);
  a[0] = 3.; a[1] = 2.;
  a[2] = 4.; a[3] = 5.;
  myvec_t b(4);
  b[0] = 1.; b[1] = 0.;
  b[2] = 2.; b[3] = 3.;
  myvec_t c(4);
  c[0] = 2.; c[1] = 1.;
  c[2] = 1.; c[3] = 3.;

  a += b - c;//a = a + (b - c)
  EXPECT_DOUBLE_EQ( a[0], 2.0);
  EXPECT_DOUBLE_EQ( a[1], 1.0);
  EXPECT_DOUBLE_EQ( a[2], 5.0);
  EXPECT_DOUBLE_EQ( a[3], 5.0);
}

TEST(algebra_vector_teuchos_serial_dense_class,
     compoundAssignMinus){

  myvec_t a(4);
  a[0] = 3.; a[1] = 2.;
  a[2] = 4.; a[3] = 5.;
  myvec_t b(4);
  b[0] = 1.; b[1] = 0.;
  b[2] = 2.; b[3] = 3.;
  myvec_t c(4);
  c[0] = 2.; c[1] = 1.;
  c[2] = 1.; c[3] = 3.;

  a -= b - c; //a = a - (b -c)
  EXPECT_DOUBLE_EQ( a[0], 4.0);
  EXPECT_DOUBLE_EQ( a[1], 3.0);
  EXPECT_DOUBLE_EQ( a[2], 3.0);
  EXPECT_DOUBLE_EQ( a[3], 5.0);

  a[0] = 3.; a[1] = 2.;
  a[2] = 4.; a[3] = 5.;
  a = a - (b - c); //a = a - (b -c)
  EXPECT_DOUBLE_EQ( a[0], 4.0);
  EXPECT_DOUBLE_EQ( a[1], 3.0);
  EXPECT_DOUBLE_EQ( a[2], 3.0);
  EXPECT_DOUBLE_EQ( a[3], 5.0);
}


TEST(algebra_vector_teuchos_serial_dense_class,
     various){

  myvec_t a(4);
  a[0] = 3.; a[1] = 2.;
  a[2] = 4.; a[3] = 5.;
  myvec_t b(4);
  b[0] = 1.; b[1] = 0.;
  b[2] = 2.; b[3] = 3.;
  myvec_t c(4);
  c[0] = 2.; c[1] = 1.;
  c[2] = 1.; c[3] = 3.;

  constexpr double f1 = 2.;
  constexpr double f2 = 3.;

  myvec_t d(4);
  d = (a - c) - b;
  EXPECT_DOUBLE_EQ( d[0], 0.0);
  EXPECT_DOUBLE_EQ( d[1], 1.0);
  EXPECT_DOUBLE_EQ( d[2], 1.0);
  EXPECT_DOUBLE_EQ( d[3], -1.0);
  d = (a - b) + c;
  EXPECT_DOUBLE_EQ( d[0], 4.0);
  EXPECT_DOUBLE_EQ( d[1], 3.0);
  EXPECT_DOUBLE_EQ( d[2], 3.0);
  EXPECT_DOUBLE_EQ( d[3], 5.0);

  d = a*f1 + c * b * f2;
  EXPECT_DOUBLE_EQ( d[0], 12.0);
  EXPECT_DOUBLE_EQ( d[1], 4.0);
  EXPECT_DOUBLE_EQ( d[2], 14.0);
  EXPECT_DOUBLE_EQ( d[3], 37.0);
  d = a*f1 + c * f2 * b;
  EXPECT_DOUBLE_EQ( d[0], 12.0);
  EXPECT_DOUBLE_EQ( d[1], 4.0);
  EXPECT_DOUBLE_EQ( d[2], 14.0);
  EXPECT_DOUBLE_EQ( d[3], 37.0);
  d = a*f1 + c * b * 3.;
  EXPECT_DOUBLE_EQ( d[0], 12.0);
  EXPECT_DOUBLE_EQ( d[1], 4.0);
  EXPECT_DOUBLE_EQ( d[2], 14.0);
  EXPECT_DOUBLE_EQ( d[3], 37.0);
  d = a*f1 + f2 * c * b;
  EXPECT_DOUBLE_EQ( d[0], 12.0);
  EXPECT_DOUBLE_EQ( d[1], 4.0);
  EXPECT_DOUBLE_EQ( d[2], 14.0);
  EXPECT_DOUBLE_EQ( d[3], 37.0);

  d = a + b + c;
  EXPECT_DOUBLE_EQ( d[0], 6.0);
  EXPECT_DOUBLE_EQ( d[1], 3.0);
  EXPECT_DOUBLE_EQ( d[2], 7.0);
  EXPECT_DOUBLE_EQ( d[3], 11.0);
  d = a + c + b;
  EXPECT_DOUBLE_EQ( d[0], 6.0);
  EXPECT_DOUBLE_EQ( d[1], 3.0);
  EXPECT_DOUBLE_EQ( d[2], 7.0);
  EXPECT_DOUBLE_EQ( d[3], 11.0);

  d = a + b * c;
  EXPECT_DOUBLE_EQ( d[0], 5.0);
  EXPECT_DOUBLE_EQ( d[1], 2.0);
  EXPECT_DOUBLE_EQ( d[2], 6.0);
  EXPECT_DOUBLE_EQ( d[3], 14.0);
  d = a + c * b;
  EXPECT_DOUBLE_EQ( d[0], 5.0);
  EXPECT_DOUBLE_EQ( d[1], 2.0);
  EXPECT_DOUBLE_EQ( d[2], 6.0);
  EXPECT_DOUBLE_EQ( d[3], 14.0);

  d = f1*a + b * c;
  EXPECT_DOUBLE_EQ( d[0], 8.0);
  EXPECT_DOUBLE_EQ( d[1], 4.0);
  EXPECT_DOUBLE_EQ( d[2], 10.0);
  EXPECT_DOUBLE_EQ( d[3], 19.0);
  d = a*f1 + c * b;
  EXPECT_DOUBLE_EQ( d[0], 8.0);
  EXPECT_DOUBLE_EQ( d[1], 4.0);
  EXPECT_DOUBLE_EQ( d[2], 10.0);
  EXPECT_DOUBLE_EQ( d[3], 19.0);
}
