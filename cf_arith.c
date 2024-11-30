//--------------------------------------------------------------------//
// continued fractions arithmetic with Gosper's algorithm (cf HAKMEM) //
//--------------------------------------------------------------------//
// the multiple precision GNU library GMP is required (cf http://www.swox.com/gmp/)

#define READLINE_LIBRARY

#include <errno.h>
#include <stdio.h> // for printf()
#include <stdlib.h> // for malloc()
#include <string>
#include <iostream>
#include <setjmp.h>
#include <signal.h>
#include <unistd.h>
#include <readline/readline.h>
#include <readline/history.h>

#include "cf_arith.h" // public interfaces of classes

//#define DEBUG
//#define YDEBUG
//#define ZDEBUG
//#define XDEBUG

//-------------------------//
// constants : 0, 1 and 10 //
//-------------------------//
mpz_t zero, one, ten;
mpz_t scale; // length of cf output

//----------------//
//----------------//
//  class Vector  //
//----------------//
//----------------//

// constructor
//--------------//
Vector::Vector() {
//--------------//
#ifdef DEBUG
   printf("Vector::Vector() %d\n", this);
#endif
   mpz_init(a); mpz_init(b);
}

// constructor
//--------------------------------//
Vector::Vector(mpz_t xa, mpz_t xb) {
//--------------------------------//
#ifdef DEBUG
   printf("Vector::Vector(mpz_t, mpz_t) %d\n", this);
#endif
   mpz_init(a); mpz_init(b);
   mpz_set(a, xa); mpz_set(b, xb);
}

// constructor
//----------------------------//
Vector::Vector(int xa, int xb) {
//----------------------------//
#ifdef DEBUG
   printf("Vector::Vector(int, int) %d\n", this);
#endif
   mpz_init(a); mpz_init(b);
   mpz_set_si(a, xa); mpz_set_si(b, xb);
}

// destructor
//---------------//
Vector::~Vector() {
//---------------//
#ifdef DEBUG
   printf("Vector::~Vector() %d\n", this);
#endif
   mpz_clear(a); mpz_clear(b);
}

// initialisation
//------------------------------//
Vector::Vector(const Vector & v) {
//------------------------------//
#ifdef DEBUG
   printf("Vector::Vector(const Vector & v) %d\n", this);
#endif
   mpz_init(a); mpz_init(b);
   mpz_set(a, v.a);
   mpz_set(b, v.b);
}

// assignment
//--------------------------------------------//
Vector & Vector::operator=(const Vector & rsh) {
//--------------------------------------------//
#ifdef DEBUG
   printf("Vector & Vector::operator=(const Vector &) %d %d\n", this, &rsh);
#endif
   if (this != &rsh) {
      mpz_set(a, rsh.a);
      mpz_set(b, rsh.b);
   }
   return *this;
}

// in place transformation
//---------------------------//
void Vector::VectorSimplify() {
//---------------------------//
   mpz_t tmp; mpz_init(tmp);
   mpz_gcd(tmp, a, b);
   if (mpz_cmp_ui(tmp, 1) > 0) {
      mpz_fdiv_q(a, a, tmp); // a /= gcd(a, b);
      mpz_fdiv_q(b, b, tmp); // b /= gcd(a, b);
   }
   mpz_clear(tmp);
}

//----------------//
//----------------//
//  class Matrix  //
//----------------//
//----------------//

// private initialization
//-----------------//
void Matrix::init() {
//-----------------//
#ifdef DEBUG
   printf("Matrix::init() %d\n", this);
#endif
   mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);
}

// constructor
//--------------//
Matrix::Matrix() {
//--------------//
#ifdef DEBUG
   printf("Matrix::Matrix() %d\n", this);
#endif
   init();
}

// constructor
//----------------------------------------------------//
Matrix::Matrix(mpz_t xa, mpz_t xb, mpz_t xc, mpz_t xd) {
//----------------------------------------------------//
#ifdef DEBUG
   printf("Matrix::Matrix(mpz_t, mpz_t, mpz_t, mpz_t) %d\n", this);
#endif
   init();
   mpz_set(a, xa); mpz_set(b, xb); mpz_set(c, xc); mpz_set(d, xd);
}

// constructor
//--------------------------------------------//
Matrix::Matrix(int xa, int xb, int xc, int xd) {
//--------------------------------------------//
#ifdef DEBUG
   printf("Matrix::Matrix(int, int, int, int) %d\n", this);
#endif
   init();
   mpz_set_si(a, xa); mpz_set_si(b, xb); mpz_set_si(c, xc); mpz_set_si(d, xd);
}

//        /     \
// m  <-  |v1 v2|
//        \     /
// constructor
//--------------------------------------//
Matrix::Matrix(Vector & v1, Vector  &v2) {
//--------------------------------------//
   init();
   mpz_set(a, *v1._geta()); mpz_set(b, *v2._geta());
   mpz_set(c, *v1._getb()); mpz_set(d, *v2._getb());
}

// destructor
//---------------//
Matrix::~Matrix() {
//---------------//
#ifdef DEBUG
   printf("Matrix::~Matrix() %d\n", this);
#endif
   mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);
}

// initialisation
//------------------------------//
Matrix::Matrix(const Matrix & m) {
//------------------------------//
#ifdef DEBUG
   printf("Matrix::Matrix(const Matrix & m) %d\n", this);
#endif
   init();
   mpz_set(a, m.a); mpz_set(b, m.b);
   mpz_set(c, m.c); mpz_set(d, m.d);
}

// assignment
//--------------------------------------------------//
const Matrix & Matrix::operator=(const Matrix & rsh) {
//--------------------------------------------------//
#ifdef DEBUG
   printf("Matrix & Matrix::operator=(const Matrix &) %d %d\n", this, &rsh);
#endif
   if (this != &rsh) {
      mpz_set(a, rsh.a); mpz_set(b, rsh.b);
      mpz_set(c, rsh.c); mpz_set(d, rsh.d);
   }
   return *this;
}

// m  +=  m1
//--------------------------------------------------//
const Matrix & Matrix::operator+=(const Matrix & m1) {
//--------------------------------------------------//
   mpz_add(a, a, m1.a);
   mpz_add(b, b, m1.b);
   mpz_add(c, c, m1.c);
   mpz_add(d, d, m1.d);
   return *this;
}

// m  =  m1 + m2
//----------------------------------------------------------//
const Matrix operator+(const Matrix & m1, const Matrix & m2) {
//----------------------------------------------------------//
   Matrix m = m1;
   return m += m2;
}

// m  *=  m1
//--------------------------------------------------//
const Matrix & Matrix::operator*=(const Matrix & m1) {
//--------------------------------------------------//
   mpz_t xa, xb, xc, xd, tmp;
   mpz_init(xa); mpz_init(xb); mpz_init(xc); mpz_init(xd); mpz_init(tmp);

   // a = a*a1 + b*c1
   mpz_mul(xa, a, m1.a);
   mpz_mul(tmp, b, m1.c);
   mpz_add(xa, xa, tmp);
   // b = a*b1 + b*d1
   mpz_mul(xb, a, m1.b);
   mpz_mul(tmp, b, m1.d);
   mpz_add(xb, xb, tmp);
   // c = c*a1 + d*c1
   mpz_mul(xc, c, m1.a);
   mpz_mul(tmp, d, m1.c);
   mpz_add(xc, xc, tmp);
   // d = c*b1 + d*d1
   mpz_mul(xd, c, m1.b);
   mpz_mul(tmp, d, m1.d);
   mpz_add(xd, xd, tmp);
   mpz_set(a, xa); mpz_set(b, xb); mpz_set(c, xc); mpz_set(d, xd);
   mpz_clear(xa); mpz_clear(xb); mpz_clear(xc); mpz_clear(xd); mpz_clear(tmp);
   return *this;
}

// m  =  m1 . m2
//----------------------------------------------------------//
const Matrix operator*(const Matrix & m1, const Matrix & m2) {
//----------------------------------------------------------//
   Matrix m = m1;
   return m *= m2;
}

// m  =  i * m1
//------------------------------------------------------//
const Matrix operator*(const mpz_t i, const Matrix & m1) {
//------------------------------------------------------//
   Matrix m;
   mpz_mul(m.a, i, m1.a);
   mpz_mul(m.b, i, m1.b);
   mpz_mul(m.c, i, m1.c);
   mpz_mul(m.d, i, m1.d);
   return m;
}

//----------------------------------------------------
//                   /ma mb\   /va\     /ma*va + mb*vb\
// v  =  m1 . v1  =  |     | . |  |  =  |             |
//                   \mc md/   \vb/     \mc*va + md*vb/
//----------------------------------------------------
//----------------------------------------------------------//
const Vector operator*(const Matrix & m1, const Vector & v1) {
//----------------------------------------------------------//
   Vector v;
   mpz_t tmp; mpz_init(tmp);

   mpz_mul(tmp, m1.a, v1.a);
   mpz_mul(v.a, m1.b, v1.b);
   mpz_add(v.a, v.a, tmp);
   mpz_mul(tmp, m1.c, v1.a);
   mpz_mul(v.b, m1.d, v1.b);
   mpz_add(v.b, v.b, tmp);
   mpz_clear(tmp);
   return v;
}

//------------------------//
void Matrix::matrixPrint() { // print (a b c d)
//------------------------//
   printf("("); mpz_out_str(stdout, 10, a); printf(" ");
                mpz_out_str(stdout, 10, b); printf("\n");
                mpz_out_str(stdout, 10, c); printf(" ");
                mpz_out_str(stdout, 10, d); printf(")\n");
}

// in place transformation
//----------------------------//
void Matrix::matrixTranspose() {
//----------------------------//
   mpz_t tmp; mpz_init(tmp);
   mpz_set(tmp, b); mpz_set(b, c); mpz_set(c, tmp);
   mpz_clear(tmp);
}

//----------------//
//----------------//
//  class Tensor  //
//----------------//
//----------------//

// constructor
//--------------//
Tensor::Tensor() {
//--------------//
#ifdef DEBUG
   printf("Tensor::Tensor() %d\n", this);
#endif
   // not necessary because of Matrix::Matrix() !!!
   // mpz_init(*matrix_1._getx()); mpz_init(*matrix_2._gety());
}

// constructor
//--------------------------------------//
Tensor::Tensor(Matrix & m1, Matrix & m2) {
//--------------------------------------//
#ifdef DEBUG
   printf("Tensor::Tensor(Matrix &, Matrix &) %d\n", this);
#endif
   matrix_1 = m1;
   matrix_2 = m2;
}

// initialisation
//------------------------------//
Tensor::Tensor(const Tensor & t) {
//------------------------------//
#ifdef DEBUG
   printf("Tensor::Tensor(const Tensor & t) %d\n", this);
#endif
   matrix_1 = t.matrix_1;
   matrix_2 = t.matrix_2;
}

// destructor
//---------------//
Tensor::~Tensor() {
//---------------//
#ifdef DEBUG
   printf("Tensor::~Tensor() %d\n", this);
#endif
   // not necessary because of Matrix::~Matrix() !!!
   // mpz_clear(*matrix_1._getx()); mpz_clear(*matrix_2._gety());
}

// assignment
//--------------------------------------------//
Tensor & Tensor::operator=(const Tensor & rsh) {
//--------------------------------------------//
#ifdef DEBUG
   printf("Tensor & Tensor::operator=(const Tensor &) %d %d\n", this, &rsh);
#endif
   if (this != &rsh) {
      matrix_1 = rsh.matrix_1;
      matrix_2 = rsh.matrix_2;
   }
   return *this;
}

//------------------------//
void Tensor::tensorPrint() { // print (a b c d) / (e f g h)
//------------------------//
   printf("("); mpz_out_str(stdout, 10, *matrix_1._geta()); printf(" ");
                mpz_out_str(stdout, 10, *matrix_1._getb()); printf(" ");
                mpz_out_str(stdout, 10, *matrix_1._getc()); printf(" ");
                mpz_out_str(stdout, 10, *matrix_1._getd()); printf(") / (");
                mpz_out_str(stdout, 10, *matrix_2._geta()); printf(" ");
                mpz_out_str(stdout, 10, *matrix_2._getb()); printf(" ");
                mpz_out_str(stdout, 10, *matrix_2._getc()); printf(" ");
                mpz_out_str(stdout, 10, *matrix_2._getd()); printf(")\n");
}

//                 //      \   /v1\\T     /         \T
// m  =  t * v  =  ||m1, m2| * |  ||   =  |m1.v m2.v|
//                 \\      /   \v2//      \         /
//--------------------------------------------------------//
const Matrix operator*(const Tensor & t, const Vector & v) {
//--------------------------------------------------------//
#ifdef DEBUG
   printf("const Matrix operator*(const Tensor & t, const Vector & v)\n");
#endif
   Vector tmp1, tmp2;
   tmp1 = t._getm1() * v;
   tmp2 = t._getm2() * v;
   Matrix m(tmp1, tmp2);
   m.matrixTranspose(); // ATTENTION !!!
   return m;
}

//                 /a b\   /      \     /                    \
// t  =  m * t  =  |   | * |m1, m2|  =  |a*m1+b*m2, c*m1+d*m2|
//                 \c d/   \      /     \                    /
//--------------------------------------------------------//
const Tensor operator*(const Matrix & m, const Tensor & t) {
//--------------------------------------------------------//
   Matrix m1, m2;
   m1 = (m.a * t.matrix_1) + (m.b * t.matrix_2);
   m2 = (m.c * t.matrix_1) + (m.d * t.matrix_2);
   Tensor t1(m1, m2);
   return t1;
}

// in place transformation
//----------------------------//
void Tensor::tensorTranspose() {
//----------------------------//
   matrix_1.matrixTranspose();
   matrix_2.matrixTranspose();
}

//                 /      \         /          \
// t  =  t * m  =  |m1, m2| * m  =  |m1.m, m2.m|
//                 \      /         \          /
//--------------------------------------------------------//
const Tensor operator*(const Tensor & t, const Matrix & m) {
//--------------------------------------------------------//
   Matrix m1, m2;
   m1 = t.matrix_1 * m;
   m2 = t.matrix_2 * m;
   Tensor t1(m1, m2);
   return t1;
}

//---------------------------//
//---------------------------//
//  class ContinuedFraction  //
//---------------------------//
//---------------------------//

// constructor
//------------------------------------//
ContinuedFraction::ContinuedFraction() {
//------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction() %d\n", this);
#endif
   //mpz_init(*v._getx()); // not necessary because of Vector();
   //mpz_init(*xunnormNumber.m._getx()); // not necessary because of Matrix();
   //mpz_init(*xgeneralizedUnnormNumber.m._getx()); // not necessary because of Matrix();
   mpz_init(xgeneralizedUnnormNumber.n);
   mpz_init(xgeneralizedUnnormNumber.d);
   //mpz_init(*xunaryOperation.m._geta()); // not necessary because of Matrix();
   // x0 and y0 init
   mpz_init(x0); mpz_init(x1);
   mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// constructor (rational)
//------------------------------------------------------------//
ContinuedFraction::ContinuedFraction(int a, int b,
                                     enumOutput xdiscriOutput) {
//------------------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction(int, int, enumOutput) %d\n", this);
#endif
   discriNode = rational;
   discriOutput = xdiscriOutput;
   mpz_set_si(*v._geta(), a);
   mpz_set_si(*v._getb(), b);
   v.VectorSimplify(); // optional
   // x0 and y0 init
   mpz_init(x0); mpz_init(x1);
   mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// constructor (rational)
//------------------------------------------------------------//
ContinuedFraction::ContinuedFraction(mpz_t a, mpz_t b,
                                     enumOutput xdiscriOutput) {
//------------------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction(mpz_t, mpz_t, enumOutput) %d\n", this);
#endif
   discriNode = rational;
   discriOutput = xdiscriOutput;
   mpz_set(*v._geta(), a);
   mpz_set(*v._getb(), b);
   v.VectorSimplify(); // optional
   // x0 and y0 init
   mpz_init(x0); mpz_init(x1);
   mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// constructor (normNumber)
//-------------------------------------------------------------//
ContinuedFraction::ContinuedFraction(normNumberFunction function,
                                     enumOutput xdiscriOutput) {
//-------------------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction(normNumberFunction, enumOutput) %d\n", this);
#endif
   discriNode = normNumber;
   discriOutput = xdiscriOutput;
   xnormNumber.function = function;
   mpz_init(xnormNumber.state);
   mpz_set_ui(xnormNumber.state, 0);
   // x0 and y0 init
   mpz_init(x0); mpz_init(x1);
   mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// constructor (normNumber2)
//------------------------------------------------------------//
ContinuedFraction::ContinuedFraction(struct chain *xpchain,
                                     enumOutput xdiscriOutput) {
//------------------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction(struct chain *, enumOutput) %d\n", this);
#endif
   discriNode = normNumber2;
   discriOutput = xdiscriOutput;
   xnormNumber2.pchain = xpchain;
   // x0 and y0 init
   mpz_init(x0); mpz_init(x1);
   mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// constructor (unnormNumber)
//---------------------------------------------------------------//
ContinuedFraction::ContinuedFraction(unnormNumberFunction function,
                                     enumOutput xdiscriOutput) {
//---------------------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction(unnormNumberFunction, enumOutput) %d\n", this);
#endif
   discriNode = unnormNumber;
   discriOutput = xdiscriOutput;
   xunnormNumber.function = function;
   // x  =  (1.x + 0)/(0.x + 1)
   Matrix x(1, 0, 0, 1);
   xunnormNumber.m = x;
   mpz_init(xunnormNumber.state);
   mpz_set_ui(xunnormNumber.state, 0);
   // x0 and y0 init
   mpz_init(x0); mpz_init(x1);
   mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// constructor (generalizedUnnormNumber)
//--------------------------------------------------------------------------//
ContinuedFraction::ContinuedFraction(generalizedUnnormNumberFunction function,
                                     mpz_t n, mpz_t d,
                                     enumOutput xdiscriOutput) {
//--------------------------------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction(generalizedUnnormNumberFunction, mpz_t, mpz_t, enumOutput) %d\n", this);
#endif
   discriNode = generalizedUnnormNumber;
   discriOutput = xdiscriOutput;
   xgeneralizedUnnormNumber.function = function;
   // x  =  (1.x + 0)/(0.x + 1)
   Matrix x(1, 0, 0, 1);
   xgeneralizedUnnormNumber.m = x;
   mpz_init(xgeneralizedUnnormNumber.n);
   mpz_init(xgeneralizedUnnormNumber.d);
   mpz_set(xgeneralizedUnnormNumber.n, n);
   mpz_set(xgeneralizedUnnormNumber.d, d);
   mpz_init(xgeneralizedUnnormNumber.state);
   mpz_set_ui(xgeneralizedUnnormNumber.state, 0);
   // x0 and y0 init
   mpz_init(x0); mpz_init(x1);
   mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// constructor (unaryOperation)
//------------------------------------------------------------//
ContinuedFraction::ContinuedFraction(Matrix & m,
                                     ContinuedFraction & cf,
                                     enumOutput xdiscriOutput) {
//------------------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction(Matrix &, ContinuedFraction &, enumOutput) %d\n", this);
#endif
   discriNode = unaryOperation;
   discriOutput = xdiscriOutput;
   xunaryOperation.m = m;
   xunaryOperation.unary = &cf;
   // x0 and y0 init
   mpz_init(x0); mpz_init(x1);
   mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// constructor (glouton = rationalOutput)
//----------------------------------------------------------//
ContinuedFraction::ContinuedFraction(Matrix & m,
                                     ContinuedFraction & cf) {
//----------------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction(Matrix &, ContinuedFraction &) %d\n", this);
#endif
   discriNode = unaryOperation;
   discriOutput = rationalOutput;
   xunaryOperation.m = m;
   xunaryOperation.unary = &cf;
   //// x0 and y0 init
   //mpz_init(x0); mpz_init(x1);
   //mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// constructor (operation)
//------------------------------------------------------------//
ContinuedFraction::ContinuedFraction(Tensor & t,
                                     ContinuedFraction & left,
                                     ContinuedFraction & right,
                                     enumOutput xdiscriOutput) {
//------------------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::ContinuedFraction(Tensor &, ContinuedFraction &, ContinuedFraction &, enumOutput) %d\n", this);
#endif
   discriNode = operation;
   discriOutput = xdiscriOutput;
   xoperation.t = t;
   xoperation.left = &left;
   xoperation.right = &right;
   // x0 and y0 init
   mpz_init(x0); mpz_init(x1);
   mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
}

// destructor
//-------------------------------------//
ContinuedFraction::~ContinuedFraction() {
//-------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::~ContinuedFraction() %d\n", this);
#endif
}

//// assignment
////-----------------------------------------------------------------------------//
//ContinuedFraction & ContinuedFraction::operator=(const ContinuedFraction & rsh) {
////-----------------------------------------------------------------------------//
//#ifdef DEBUG
//   printf("ContinuedFraction & ContinuedFraction::operator=(const ContinuedFraction &) %d %d\n", this, &rsh);
//#endif
//   if (this != &rsh) {
//   }
//   return *this;
//}

// initialization
//-------------------------------------//
void ContinuedFraction::x0Andx1Init() {
//-------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::x0Andx1Init() %d\n", this);
#endif
   // les deux premières valeurs retournées sont poubellisées
   (void *)eval(); (void *)eval();
}

//-----------------------------------------------------------------//
ContinuedFraction ContinuedFraction::plus(ContinuedFraction & left,
                                          ContinuedFraction & right,
                                          enumOutput xdiscriOutput) {
//-----------------------------------------------------------------//
   // x + y  =  (0xy + 1x + 1y + 0)/(0xy + 0x + 0y + 1)
   Matrix m1(0, 1, 1, 0), m2(0, 0, 0, 1);
   Tensor t(m1, m2);
   return ContinuedFraction(t, left, right, xdiscriOutput);
}

//----------------------------------------------------------------------//
ContinuedFraction ContinuedFraction::substract(ContinuedFraction & left,
                                               ContinuedFraction & right,
                                               enumOutput xdiscriOutput) {
//----------------------------------------------------------------------//
   // x - y  =  (0xy + 1x - 1y + 0)/(0xy + 0x + 0y + 1)
   Matrix m1(0, 1, -1, 0), m2(0, 0, 0, 1);
   Tensor t(m1, m2);
   return ContinuedFraction(t, left, right, xdiscriOutput);
}

//---------------------------------------------------------------------//
ContinuedFraction ContinuedFraction::multiply(ContinuedFraction & left,
                                              ContinuedFraction & right,
                                              enumOutput xdiscriOutput) {
//---------------------------------------------------------------------//
   // x * y  =  (1xy + 0x + 0y + 0)/(0xy + 0x + 0y + 1)
   Matrix m1(1, 0, 0, 0), m2(0, 0, 0, 1);
   Tensor t(m1, m2);
   return ContinuedFraction(t, left, right, xdiscriOutput);
}

//-------------------------------------------------------------------//
ContinuedFraction ContinuedFraction::divide(ContinuedFraction & left,
                                            ContinuedFraction & right,
                                            enumOutput xdiscriOutput) {
//-------------------------------------------------------------------//
   // x / y  =  (0xy + 1x + 0y + 0)/(0xy + 0x + 1y + 0)
   Matrix m1(0, 1, 0, 0), m2(0, 0, 1, 0);
   Tensor t(m1, m2);
   return ContinuedFraction(t, left, right, xdiscriOutput);
}

// sqrt(b/c) initialisation
//--------------------------------------------------------------------------//
ContinuedFraction ContinuedFraction::squareRootRat(mpz_t b, mpz_t c,
                                                   enumOutput xdiscriOutput) {
//--------------------------------------------------------------------------//
#ifdef ZDEBUG
   printf("mpz_t *ContinuedFraction::squareRootRat(b,c) %d\n", this);
#endif
   discriNode = squareRoot;
   discriOutput = xdiscriOutput;
   //      /0 b\
   // m <- |   |
   //      \c 0/
   xsquareRoot.m = Matrix(zero, b, c, zero);
   // evaluation of floor(sqrt(b*c))
   mpz_t old_x; mpz_init(old_x);
   mpz_t x; mpz_init(x);
   mpz_t new_x; mpz_init(new_x);
   mpz_t det; mpz_init(det);
   mpz_t tmp1; mpz_init(tmp1);
   mpz_t tmp2; mpz_init(tmp2);
   mpz_set_ui(old_x, 0);
   mpz_set_ui(x, 0);
   mpz_set_ui(new_x, 1);
   // det = b*c
   mpz_mul(det, b, c);
   // small newton iteration : new x <- (x + det/x) / 2
   while (1) {
#ifdef ZDEBUG
      mpz_out_str(stdout, 10, new_x);
      printf(" ");
#endif
      // case x == new_x
      if (mpz_cmp(x, new_x) == 0) {
         break; // exit from while (1)
      }
      // case x == old_x+1 and new_x == old_x
      if ((mpz_cmp(old_x, new_x) == 0) &&
          (mpz_cmp(x, new_x) > 0)) {
         break; // exit from while (1)
      }
      mpz_set(old_x, x); // shift
      mpz_set(x, new_x); // shift
      mpz_tdiv_q(tmp1, det, x); // truncate (tdiv) better than floor (fdiv) ???
      mpz_add(tmp2, x, tmp1);
      mpz_tdiv_q_ui(new_x, tmp2, 2); // truncate (tdiv) better than floor (fdiv) ???
   }
#ifdef ZDEBUG
   printf("\n");
#endif
   // sqrt_det = floor(sqrt(b*c)) = floor(sqrt(det))
   mpz_init(xsquareRoot.sqrt_det);
   mpz_set(xsquareRoot.sqrt_det, new_x);
#ifdef ZDEBUG
   printf("det = ");
   mpz_out_str(stdout, 10, new_x);
   printf("\n");
#endif
   //--------------------------------------
   // periodicity detection : init the list
   //--------------------------------------
   xsquareRoot.pcurrent = NULL;
   xsquareRoot.pbegin = NULL;
   mpz_clear(old_x); mpz_clear(x); mpz_clear(new_x);
   mpz_clear(det); mpz_clear(tmp1); mpz_clear(tmp2);
   // not necessary because squareRootRat() isn't a constructor
   //// x0 and y0 init
   //mpz_init(x0); mpz_init(x1);
   //mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
   return *this;
}

//--------------------------------------------------------------------------//
ContinuedFraction ContinuedFraction::squareRootRat(mpz_t a, mpz_t b, mpz_t c,
                                                   enumOutput xdiscriOutput) {
//--------------------------------------------------------------------------//
   discriNode = squareRoot;
   discriOutput = xdiscriOutput;
   //      /a  b\
   // m <- |    |
   //      \c -a/
   mpz_t ma; mpz_init(ma); mpz_neg(ma, a);
   xsquareRoot.m = Matrix(a, b, c, ma);
   mpz_clear(ma);
   // evaluation of floor(sqrt(a*a + b*c)) = floor(sqrt(det))
   mpz_t old_x; mpz_init(old_x);
   mpz_t x; mpz_init(x);
   mpz_t new_x; mpz_init(new_x);
   mpz_t det; mpz_init(det);
   mpz_t tmp1; mpz_init(tmp1);
   mpz_t tmp2; mpz_init(tmp2);
   mpz_set_ui(old_x, 0);
   mpz_set_ui(x, 0);
   mpz_set_ui(new_x, 1);
   mpz_mul(tmp1, a, a);
   mpz_mul(tmp2, b, c);
   // det = a*a + b*c
   mpz_add(det, tmp1, tmp2);
   // small newton iteration : new x <- (x + det/x) / 2
   while (1) {
#ifdef DEBUG
      mpz_out_str(stdout, 10, new_x);
      printf(" ");
#endif
      // case x == new_x
      if (mpz_cmp(x, new_x) == 0) {
         break; // exit from while (1)
      }
      // case x == old_x+1 and new_x == old_x
      if ((mpz_cmp(old_x, new_x) == 0) &&
          (mpz_cmp(x, new_x) > 0)) {
         break; // exit from while (1)
      }
      mpz_set(old_x, x); // shift
      mpz_set(x, new_x); // shift
      mpz_tdiv_q(tmp1, det, x); // truncate (tdiv) better than floor (fdiv) ???
      mpz_add(tmp2, x, tmp1);
      mpz_tdiv_q_ui(new_x, tmp2, 2); // truncate (tdiv) better than floor (fdiv) ???
   }
#ifdef DEBUG
   printf("\n");
#endif
   // sqrt_det = floor(sqrt(a*a + b*c)) = floor(sqrt(det))
   mpz_init(xsquareRoot.sqrt_det);
   mpz_set(xsquareRoot.sqrt_det, new_x);
   //--------------------------------------
   // periodicity detection : init the list
   //--------------------------------------
   xsquareRoot.pcurrent = NULL;
   xsquareRoot.pbegin = NULL;
   mpz_clear(old_x); mpz_clear(x); mpz_clear(new_x);
   mpz_clear(det); mpz_clear(tmp1); mpz_clear(tmp2);
   // not necessary because squareRootRat() isn't a constructor
   //// x0 and y0 init
   //mpz_init(x0); mpz_init(x1);
   //mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
   return *this;
}

//-------------------------------------------------------------------------//
ContinuedFraction ContinuedFraction::squareRootOp(ContinuedFraction & cf,
                                                  enumOutput xdiscriOutput) {
//-------------------------------------------------------------------------//
   discriNode = squareRootOperation;
   discriOutput = xdiscriOutput;
   // if x = sqrt(p) and y = p
   // then x and y verifies :
   //
   //               y    0xy + 0x + 1y + 0
   //        x  =  --- = -----------------
   //               x    0xy + 1x + 0y + 0
   //
   Matrix m1(0, 0, 1, 0), m2(0, 1, 0, 0);
   Tensor t(m1, m2);
   xsquareRootOperation.t = t;
   xsquareRootOperation.unary = &cf;
   // not necessary because squareRootOp() isn't a constructor
   //// x0 and y0 init
   //mpz_init(x0); mpz_init(x1);
   //mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
   return *this;
}

// coth(1/2)
//------------------------------------------------------------------------//
ContinuedFraction ContinuedFraction::cothHalfOne(ContinuedFraction & cf,
                                                 enumOutput xdiscriOutput) {
//------------------------------------------------------------------------//
   discriNode = squareRootOperation;
   discriOutput = xdiscriOutput;
   // coth(1/2) evaluation, knowing coth(1)
   // if x = coth(1/2) and y = coth(1)
   // then x and y verifies :
   //
   //              x y - 1   1xy + 0x + 0y - 1
   //        x  =  ------- = -----------------
   //               x - y    0xy + 1x - 1y + 0
   //
   Matrix m1(1, 0, 0, -1), m2(0, 1, -1, 0);
   Tensor t(m1, m2);
   xsquareRootOperation.t = t;
   xsquareRootOperation.unary = &cf;
   // not necessary because cothHalfOne() isn't a constructor
   //// x0 and y0 init
   //mpz_init(x0); mpz_init(x1);
   //mpz_set_ui(x0, 0); mpz_set_ui(x1, 0);
   return *this;
}

//------------
// evaluators
//------------

// generic evaluator
//------------------------------//
mpz_t *ContinuedFraction::eval() {
//------------------------------//
#ifdef DEBUG
   printf("mpz_t *ContinuedFraction::eval() %d\n", this);
#endif
   switch (discriNode) {
      case rational:
         return eval_rational();
      case normNumber:
         return eval_normNumber();
      case normNumber2:
         return eval_normNumber2();
      case unnormNumber:
         return eval_unnormNumber();
      case generalizedUnnormNumber:
         return eval_generalizedUnnormNumber();
      case operation:
         return eval_operation();
      case unaryOperation:
         return eval_unaryOperation();
      case squareRoot:
         return eval_squareRoot();
      case squareRootOperation:
         return eval_squareRootOperation();
      default:
         printf("mpz_t *ContinuedFraction::eval() ERROR\n");
         exit(-1);
   }
}

// rational evaluator
//---------------------------------------//
mpz_t *ContinuedFraction::eval_rational() {
//---------------------------------------//
#ifdef DEBUG
   printf("mpz_t *ContinuedFraction::eval_rational() %d\n", this);
   printf("("); mpz_out_str(stdout, 10, *v._geta()); printf(" ");
                mpz_out_str(stdout, 10, *v._getb()); printf(")\n");
#endif
   mpz_t *pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result
   mpz_t quo1, tmp; mpz_init(quo1); mpz_init(tmp);

   mpz_set(*pquo1, x0); // shift
   mpz_set(x0, x1); // shift
   if (mpz_cmp_ui(*v._getb(), 0) == 0) { // end of rational
#ifdef YDEBUG
   printf("\nrational -> end\n");
#endif
      // ce n'est pas genant de passer trois fois dans ce cas a la fin du developpement
      //mpz_set_ui(x1, 0); // end of rational
      mpz_set_si(x1, -77); // end of rational
   } else {
      mpz_fdiv_q(quo1, *v._geta(), *v._getb()); // floor (fdiv) 
      //-------------
      // output quo1
      //-------------
      mpz_neg(quo1, quo1);
      if (discriOutput == cfOutput) {
#if 1
         Matrix moutput(zero, one, one, quo1);
         v = moutput * v;
#else
         // the two lines below are equivalent to m * v
         mpz_set(v.a, v.b);
         mpz_set(v.b, tmp);
#endif
      } else if (discriOutput == decimalOutput) {
         mpz_mul_ui(tmp, quo1, 10);
         Matrix moutput(ten, tmp, zero, one);
         v = moutput * v;
      }
      mpz_neg(quo1, quo1);
      mpz_set(x1, quo1); // store result
#ifdef DEBUG
      printf("ContinuedFraction::eval_rational() quo1 = ");
      mpz_out_str(stdout, 10, quo1); printf("\n");
      printf("ContinuedFraction::eval_rational() a = ");
      mpz_out_str(stdout, 10, *v._geta()); printf("\n");
      printf("ContinuedFraction::eval_rational() b = ");
      mpz_out_str(stdout, 10, *v._getb()); printf("\n");
#endif
   }
   mpz_clear(quo1); mpz_clear(tmp);
   return pquo1; // that's all folks
}

// normNumber evaluator
//-----------------------------------------//
mpz_t *ContinuedFraction::eval_normNumber() {
//-----------------------------------------//
#ifdef DEBUG
   printf("mpz_t *ContinuedFraction::eval_normNumber() %d\n", this);
#endif
   mpz_t *pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result
   mpz_t *eval;

   //-------------
   // output quo1
   //-------------
   if (discriOutput == cfOutput) {
      mpz_set(*pquo1, x0);
      mpz_set(x0, x1);
      //mpz_set(x1, *((*xnormNumber.function)(xnormNumber.state++)));
      eval = (*xnormNumber.function)(xnormNumber.state);
      mpz_set(x1, *eval);
      mpz_clear(*eval);
      mpz_add_ui(xnormNumber.state, xnormNumber.state, 1);
#ifdef DEBUG
      mpz_out_str(stdout, 10, *pquo1); printf("\n");
#endif
      return pquo1; // that's all folks
   } else if (discriOutput == decimalOutput) {
      printf("impossible error...\n");
      exit(-1);
   }
}

// normNumber2 evaluator
//------------------------------------------//
mpz_t *ContinuedFraction::eval_normNumber2() {
//------------------------------------------//
#ifdef DEBUG
   printf("mpz_t *ContinuedFraction::eval_normNumber2() %d\n", this);
#endif
   mpz_t *pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result

   //-------------
   // output quo1
   //-------------
   if (discriOutput == cfOutput) {
      // shift
      mpz_set(*pquo1, x0);
      mpz_set(x0, x1);
      if (xnormNumber2.pchain != NULL) {
         mpz_set(x1, xnormNumber2.pchain->value);
         xnormNumber2.pchain = xnormNumber2.pchain->next;
      } else {
         // end of normNumber2
         //mpz_set_ui(x1, 0);
         mpz_set_si(x1, -77);
      }
#ifdef DEBUG
      mpz_out_str(stdout, 10, *pquo1); printf("\n");
#endif
      return pquo1; // that's all folks
   } else if (discriOutput == decimalOutput) {
      printf("impossible error...\n");
      exit(-1);
   }
}

// unnormNumber evaluator
//-------------------------------------------//
mpz_t *ContinuedFraction::eval_unnormNumber() {
//-------------------------------------------//
#ifdef DEBUG
   printf("mpz_t *ContinuedFraction::eval_unnormNumber() %d\n", this);
   xunnormNumber.m.matrixPrint();
#endif
   int xoutput;
   Vector xvector;
   mpz_t quo1, quo2, tmp;
   mpz_t *pquo1;

   // get next vector input
   //xvector = (*xunnormNumber.function)(xunnormNumber.state++);
   xvector = (*xunnormNumber.function)(xunnormNumber.state);
   mpz_add_ui(xunnormNumber.state, xunnormNumber.state, 1);
   if (mpz_cmp_ui(*xvector._getb(), 0) == 0) { // end of unnormNumber
      //--------------------------
      // unnormNumber -> rational
      //--------------------------
#ifdef YDEBUG
      printf("\nunnormnumber -> rational\n");
#endif
      discriNode = rational;
      Vector v2(*xvector._geta(), one);
      v = xunnormNumber.m * v2;
      v.VectorSimplify();
#ifdef DEBUG
      printf("unnormNumber -> rational\n");
      printf(" vector ("); mpz_out_str(stdout, 10, *xvector._geta()); printf(" ");
                           mpz_out_str(stdout, 10, one); printf(")\n");
      printf("("); mpz_out_str(stdout, 10, *v._geta()); printf(" ");
                   mpz_out_str(stdout, 10, *v._getb()); printf(")\n");
#endif
      return eval(); // recursive call
   } else {
      mpz_init(quo1); mpz_init(quo2); mpz_init(tmp);
      if ((mpz_cmp_ui(*xunnormNumber.m._getc(), 0) != 0) &&
          (mpz_cmp_ui(*xunnormNumber.m._getd(), 0) != 0)) {
         mpz_fdiv_q(quo1, *xunnormNumber.m._geta(), // floor (fdiv)
                          *xunnormNumber.m._getc());
         mpz_fdiv_q(quo2, *xunnormNumber.m._getb(), // floor (fdiv)
                          *xunnormNumber.m._getd());
         if (mpz_cmp(quo1, quo2) == 0) {
            xoutput = 1; // output
         } else {
            xoutput = 0; // input
         }
      } else {
         xoutput = 0; // input
      }
      if (xoutput) {
         // unget next vector input
         mpz_sub_ui(xunnormNumber.state, xunnormNumber.state, 1);
         //-------------
         // output quo1
         //-------------
         mpz_neg(quo1, quo1);
         if (discriOutput == cfOutput) {
            Matrix moutput(zero, one, one, quo1);
            xunnormNumber.m = moutput * xunnormNumber.m;
         } else if (discriOutput == decimalOutput) {
            mpz_mul_ui(tmp, quo1, 10);
            Matrix moutput(ten, tmp, zero, one);
            xunnormNumber.m = moutput * xunnormNumber.m;
         }
         mpz_neg(quo1, quo1);
#ifdef DEBUG
         printf("unnormNumber : output ");
         mpz_out_str(stdout, 10, quo1); printf("\n");
         xunnormNumber.m.matrixPrint();
#endif
         pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result
         mpz_set(*pquo1, x0);
         mpz_set(x0, x1);
         mpz_set(x1, quo1);
         mpz_clear(quo1); mpz_clear(quo2); mpz_clear(tmp);
         return pquo1; // that's all folks
      } else {
         //------------
         // input a, b
         //------------
         Matrix minput(*xvector._geta(), *xvector._getb(), one, zero);
         xunnormNumber.m = xunnormNumber.m * minput;
#ifdef DEBUG
         printf("unnormNumber : input (");
         mpz_out_str(stdout, 10, *xvector._geta()); printf(" ");
         mpz_out_str(stdout, 10, *xvector._getb()); printf(")\n");
         xunnormNumber.m.matrixPrint();
#endif
         mpz_clear(quo1); mpz_clear(quo2); mpz_clear(tmp);
         return eval(); // recursive call
      }
   }
}

// generalizedUnnormNumber evaluator
//------------------------------------------------------//
mpz_t *ContinuedFraction::eval_generalizedUnnormNumber() {
//------------------------------------------------------//
#ifdef XDEBUG
   printf("mpz_t *ContinuedFraction::eval_generalizedUnnormNumber() %d\n", this);
   //xgeneralizedUnnormNumber.m.matrixPrint();
#endif
   int xoutput;
   Matrix xmatrix;
   mpz_t *pquo1;
   mpz_t quo1, quo2, quo3, tmp1, tmp2;

   // get next matrix input
   xmatrix = (*xgeneralizedUnnormNumber.function)(xgeneralizedUnnormNumber.state,
                                                  xgeneralizedUnnormNumber.n,
                                                  xgeneralizedUnnormNumber.d);
   mpz_add_ui(xgeneralizedUnnormNumber.state, xgeneralizedUnnormNumber.state, 1);
// A MODIFIER !!!
#if 0 /*{*/
   if (mpz_cmp_ui(*xmatrix._getc(), 0) == 0) { // end of generalizedUnnormNumber
      //-------------------------------------
      // generalizedUnnormNumber -> rational
      //-------------------------------------
#ifdef YDEBUG
      printf("\ngeneralizedunnormnumber -> rational\n");
#endif
      discriNode = rational;
      Vector v2(*xmatrix._geta(), one);
      v = xgeneralizedUnnormNumber.m * v2;
      v.VectorSimplify();
#ifdef DEBUG
      printf("generalizedUnnormnumber -> rational\n");
      printf("("); mpz_out_str(stdout, 10, *v._geta()); printf(" ");
                   mpz_out_str(stdout, 10, *v._getb()); printf(")\n");
#endif
      return eval(); // recursive call
   } else {
#endif /*}}{*/
      mpz_init(quo1); mpz_init(quo2); mpz_init(quo3); mpz_init(tmp1); mpz_init(tmp2);
      mpz_add(tmp1, *xgeneralizedUnnormNumber.m._geta(), *xgeneralizedUnnormNumber.m._getb());
      mpz_add(tmp2, *xgeneralizedUnnormNumber.m._getc(), *xgeneralizedUnnormNumber.m._getd());
      // on est obligé de faire les calculs alambiqués ci-dessous à cause de tan(3) et tan(178/113)...
      // sinon on peut avoir des quotients négatifs (ou nuls) en cours de développement...
      // le test floor(a/c) = floor(b/d) n'est pas suffisant...
      // le test floor(a/c) = floor((a+b)/(c+d)) n'est pas suffisant contrairement à l'article de I. Vardi...
      // le test complet est floor(a/c) = floor(b/d) ET floor(a/c) = floor((a+b)/(c+d))
      // on peut avoir des quotients nuls en cours de développement : log(1/2) ou log(1/3)...
      // (le développement n'est pas vraiment convergent pour x<1)
      // ces cas sont éliminés par la formule log(1/x) = -log(x)
      if ((mpz_cmp_ui(*xgeneralizedUnnormNumber.m._getc(), 0) != 0) && // c != 0 && d != 0 && c+d != 0
          (mpz_cmp_ui(*xgeneralizedUnnormNumber.m._getd(), 0) != 0) &&
          (mpz_cmp_ui(tmp2, 0) != 0)) {
         mpz_fdiv_q(quo1, *xgeneralizedUnnormNumber.m._geta(), // floor (fdiv)
                    *xgeneralizedUnnormNumber.m._getc());
         mpz_fdiv_q(quo2, tmp1, tmp2); // floor (fdiv)
         mpz_fdiv_q(quo3, *xgeneralizedUnnormNumber.m._getb(), // floor (fdiv)
                    *xgeneralizedUnnormNumber.m._getd());
         // on teste floor(a/c) = floor((a+b)/(c+d)) = floor(b/d)
         if ((mpz_cmp(quo1, quo2) == 0) && (mpz_cmp(quo1, quo3) == 0)) {
            xoutput = 1; // output
         } else {
            xoutput = 0; // input
         }
      } else {
         xoutput = 0; // input
      }
      if (xoutput) {
         // unget next matrix input
         mpz_sub_ui(xgeneralizedUnnormNumber.state, xgeneralizedUnnormNumber.state, 1);
         //-------------
         // output quo1
         //-------------
         mpz_neg(quo1, quo1);
         if (discriOutput == cfOutput) {
            Matrix moutput(zero, one, one, quo1);
            xgeneralizedUnnormNumber.m = moutput * xgeneralizedUnnormNumber.m;
         } else if (discriOutput == decimalOutput) {
            mpz_mul_ui(tmp1, quo1, 10);
            Matrix moutput(ten, tmp1, zero, one);
            xgeneralizedUnnormNumber.m = moutput * xgeneralizedUnnormNumber.m;
         }
         mpz_neg(quo1, quo1);
#ifdef XDEBUG
         printf("generalizedUnnormNumber : output : "); mpz_out_str(stdout, 10, quo1); printf("\n");
         xgeneralizedUnnormNumber.m.matrixPrint();
#endif
         pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result
         mpz_set(*pquo1, x0);
         mpz_set(x0, x1);
         mpz_set(x1, quo1);
         mpz_clear(quo1); mpz_clear(quo2); mpz_clear(quo3); mpz_clear(tmp1); mpz_clear(tmp2);
         return pquo1; // that's all folks
      } else {
         //-------------------------
         // input matrix a, b, c, d
         //-------------------------
         xgeneralizedUnnormNumber.m = xgeneralizedUnnormNumber.m * xmatrix;
#ifdef XDEBUG
         printf("generalizedUnnormnumber : input\n");
         xgeneralizedUnnormNumber.m.matrixPrint();
#endif
         mpz_clear(quo1); mpz_clear(quo2); mpz_clear(quo3); mpz_clear(tmp1); mpz_clear(tmp2);
         return eval(); // recursive call
      }
   //}
}

// unaryOperation evaluator
//---------------------------------------------//
mpz_t *ContinuedFraction::eval_unaryOperation() {
//---------------------------------------------//
#ifdef DEBUG
   printf("mpz_t *ContinuedFraction::eval_unaryOperation() %d\n", this);
   xunaryOperation.m.matrixPrint();
#endif
   int xoutput;
   mpz_t *pquo1;
   mpz_t quo1, quo2, tmp;

   // end of unaryOperation
   //if (mpz_cmp_ui(xunaryOperation.unary->x1, 0) == 0) {
   if (mpz_cmp_si(xunaryOperation.unary->x1, -77) == 0) {
      //----------------------------
      // unaryOperation -> rational
      //----------------------------
#ifdef YDEBUG
      printf("\nunaryoperation -> rational\n");
#endif
#ifdef DEBUG
      printf("unaryOperation -> rational ");
      mpz_out_str(stdout, 10, xunaryOperation.unary->x0);
      printf("\n");
      //printf("("); mpz_out_str(stdout, 10, *v._geta()); printf(" ");
      //             mpz_out_str(stdout, 10, *v._getb()); printf(")\n");
#endif
      discriNode = rational;
      Vector v2(xunaryOperation.unary->x0, one);
      v = xunaryOperation.m * v2;
      v.VectorSimplify();
      return eval(); // recursive call
   } else {
      mpz_init(quo1); mpz_init(quo2); mpz_init(tmp);
      if ((mpz_cmp_ui(*xunaryOperation.m._getc(), 0) != 0) &&
          (mpz_cmp_ui(*xunaryOperation.m._getd(), 0) != 0)) {
         mpz_fdiv_q(quo1, *xunaryOperation.m._geta(), // floor (fdiv)
                          *xunaryOperation.m._getc());
         mpz_fdiv_q(quo2, *xunaryOperation.m._getb(), // floor (fdiv)
                          *xunaryOperation.m._getd());
         if (mpz_cmp(quo1, quo2) == 0) {
            xoutput = 1; // output
         } else {
            xoutput = 0; // input
         }
      } else {
         xoutput = 0; // input
      }
      if (xoutput) {
         //-------------
         // output quo1
         //-------------
         mpz_neg(quo1, quo1);
         if (discriOutput == cfOutput) {
            Matrix moutput(zero, one, one, quo1);
            xunaryOperation.m = moutput * xunaryOperation.m;
         } else if (discriOutput == decimalOutput) {
            mpz_mul_ui(tmp, quo1, 10);
            Matrix moutput(ten, tmp, zero, one);
            xunaryOperation.m = moutput * xunaryOperation.m;
         } else if (discriOutput == rationalOutput) {
            printf("mpz_t *ContinuedFraction::eval_unaryOperation() ERROR\n");
            exit(-1);
         }
         mpz_neg(quo1, quo1);
#ifdef DEBUG
         printf("unaryOperation : output ");
         mpz_out_str(stdout, 10, quo1); printf("\n");
         printf("x0 = "); mpz_out_str(stdout, 10, x0); printf("\n");
         printf("x1 = "); mpz_out_str(stdout, 10, x1); printf("\n");
         //matrixPrint(xunaryOperation.m);
#endif
         pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result
         mpz_set(*pquo1, x0);
         mpz_set(x0, x1);
         mpz_set(x1, quo1);
         mpz_clear(quo1); mpz_clear(quo2); mpz_clear(tmp);
         return pquo1; // that's all folks
      } else {
         //---------
         // input a
         //---------
         Matrix minput(*(xunaryOperation.unary->eval()), one, one, zero); // recursive call
         // free(result of unary->eval); ???
         xunaryOperation.m = xunaryOperation.m * minput;
#ifdef DEBUG
         printf("unaryOperation : input ");
         //mpz_out_str(stdout, 10, quo1); printf("\n");
         //matrixPrint(xunaryOperation.m);
#endif
         mpz_clear(quo1); mpz_clear(quo2); mpz_clear(tmp);
         return eval(); // recursive call
      }
   }
}

// operation evaluator
//----------------------------------------//
mpz_t *ContinuedFraction::eval_operation() {
//----------------------------------------//
   int xoutput;
   mpz_t *pquo1;
   mpz_t quo1, quo2, quo3, quo4, tmp;

#ifdef DEBUG
   printf("ContinuedFraction::eval_operation() %d\n", this);
   xoperation.t.tensorPrint();
#endif
   // if end of right
   //if (mpz_cmp_ui(xoperation.right->x1, 0) == 0) {
   if (mpz_cmp_si(xoperation.right->x1, -77) == 0) {
      //-----------------------------
      // operation -> unaryOperation
      //-----------------------------
#ifdef YDEBUG
      printf("\noperation -> unaryoperation\n");
#endif
      discriNode = unaryOperation;
      Vector v2(xoperation.right->x0, one);
      xunaryOperation.m = xoperation.t * v2;
      xunaryOperation.unary = xoperation.left;
      return eval(); // recursive call
   } else {
      mpz_init(quo1); mpz_init(quo2); mpz_init(quo3); mpz_init(quo4);
      mpz_init(tmp);
      if ((mpz_cmp_ui(*xoperation.t._getm2()._geta(), 0) != 0) && // e
          (mpz_cmp_ui(*xoperation.t._getm2()._getb(), 0) != 0) && // f
          (mpz_cmp_ui(*xoperation.t._getm2()._getc(), 0) != 0) && // g
          (mpz_cmp_ui(*xoperation.t._getm2()._getd(), 0) != 0)) { // h
         mpz_fdiv_q(quo1, *xoperation.t._getm1()._geta(), // floor (fdiv)
                          *xoperation.t._getm2()._geta()); // e
         mpz_fdiv_q(quo2, *xoperation.t._getm1()._getb(), // floor (fdiv)
                          *xoperation.t._getm2()._getb()); // f
         mpz_fdiv_q(quo3, *xoperation.t._getm1()._getc(), // floor (fdiv)
                          *xoperation.t._getm2()._getc()); // g
         mpz_fdiv_q(quo4, *xoperation.t._getm1()._getd(), // floor (fdiv)
                          *xoperation.t._getm2()._getd()); // h
         if ((mpz_cmp(quo1, quo2) == 0) &&
             (mpz_cmp(quo2, quo3) == 0) &&
             (mpz_cmp(quo3, quo4) == 0)) {
            xoutput = 1; // output
         } else {
            xoutput = 0; // input
         }
      } else {
         xoutput = 0; // input
      }
      if (xoutput) {
#ifdef DEBUG
         printf("ContinuedFraction::eval_operation() : xoutput ");
         printf("\n"); mpz_out_str(stdout, 10, quo1); printf("\n");
         printf("x0 = "); mpz_out_str(stdout, 10, x0); printf("\n");
         printf("x1 = "); mpz_out_str(stdout, 10, x1); printf("\n");
#endif
         //-------------
         // output quo1
         //-------------
         mpz_neg(quo1, quo1);
         if (discriOutput == cfOutput) {
            Matrix moutput(zero, one, one, quo1);
            xoperation.t = moutput * xoperation.t;
         } else if (discriOutput == decimalOutput) {
            mpz_mul_ui(tmp, quo1, 10);
            Matrix moutput(ten, tmp, zero, one);
            xoperation.t = moutput * xoperation.t;
         }
         mpz_neg(quo1, quo1);
         pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result
         mpz_set(*pquo1, x0);
         mpz_set(x0, x1);
         mpz_set(x1, quo1);
         mpz_clear(quo1); mpz_clear(quo2); mpz_clear(quo3); mpz_clear(quo4);
         mpz_clear(tmp);
         return pquo1; // that's all folks
      } else {
#ifdef DEBUG
         printf("ContinuedFraction::eval_operation() : yinput and exchange\n");
#endif
         //------------------------------
         // input y and exchange x and y
         //------------------------------
         // one tensorTranspose because of exchange x and y
         Matrix m(*xoperation.right->eval(), one, one, zero); // recursive call
         // free(result of unary->eval); ???
         xoperation.t = xoperation.t * m;
         xoperation.t.tensorTranspose();
         {  ContinuedFraction *tmp;
            // left and right exchange
            tmp = xoperation.left;
            xoperation.left = xoperation.right;
            xoperation.right = tmp;
         }
         mpz_clear(quo1); mpz_clear(quo2); mpz_clear(quo3); mpz_clear(quo4);
         mpz_clear(tmp);
         return eval(); // recursive call
      }
   }
}

// squareRoot evaluator : square root of rational
//-----------------------------------------//
mpz_t *ContinuedFraction::eval_squareRoot() {
//-----------------------------------------//
   mpz_t *pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result
   mpz_t quo1, rem; mpz_init(quo1); mpz_init(rem);
#ifdef DEBUG
   printf("ContinuedFraction::eval_squareRoot() %d\n", this);
   xsquareRoot.m.matrixPrint();
   printf("det = "); mpz_out_str(stdout, 10, xsquareRoot.sqrt_det); printf("\n");
#endif
   //------------------------
   //                  /a  b\
   // m is of the form |    |
   //                  \c -a/
   //------------------------
   // if end of squareRoot (c == 0)
   if (mpz_cmp_ui(*xsquareRoot.m._getc(), 0) == 0) {
      mpz_set(*pquo1, x0);
      mpz_set(x0, x1);
      mpz_set_ui(x1, 0);
      mpz_clear(quo1); mpz_clear(rem);
#ifdef DEBUG
      printf("eval_squareRoot() end (m.c == 0) return = ");
      mpz_out_str(stdout, 10, *pquo1); printf("\n");
#endif
      return pquo1; // that's all folks
   } else {
      if (discriOutput == decimalOutput) {
         printf("impossible error...\n");
         exit(-1);
      }
      //----------------------------------------------------
      // evaluation of quo1 = floor((a + sqrt(a*a+b*c)) / c)
      //                    = floor((a + sqrt_det) / c)
      //----------------------------------------------------
      mpz_add(rem, xsquareRoot.sqrt_det, *xsquareRoot.m._geta());
      mpz_tdiv_q(quo1, rem, *xsquareRoot.m._getc());
      //--------------------------------------
      // periodicity detection : scan the list
      //--------------------------------------
      struct chain2 *pt = xsquareRoot.pbegin;
      while (pt != NULL) {
         // TODO : add the == operator on the Matrix class...
         if ((mpz_cmp(*pt->Matrixvalue._geta(), *xsquareRoot.m._geta()) == 0) &&
             (mpz_cmp(*pt->Matrixvalue._getb(), *xsquareRoot.m._getb()) == 0) &&
             (mpz_cmp(*pt->Matrixvalue._getc(), *xsquareRoot.m._getc()) == 0) &&
             (mpz_cmp(*pt->Matrixvalue._getd(), *xsquareRoot.m._getd()) == 0)) {
            // periodicity detected !!!
            // printf("\n*** periodicity detected !!! ***\n");
            // squareRoot -> normNumber2
            break; // exit from while
         }
         pt = pt->next;
      }
      //----------------------------------------------------
      // periodicity detection : append new matrix/new value
      //----------------------------------------------------
      if (xsquareRoot.pcurrent == NULL) {
         xsquareRoot.pcurrent = new struct chain2;
         xsquareRoot.pbegin = xsquareRoot.pcurrent;
      } else {
         xsquareRoot.pcurrent->next = new struct chain2;
         xsquareRoot.pcurrent = xsquareRoot.pcurrent->next;
      }
      xsquareRoot.pcurrent->next = NULL;
      xsquareRoot.pcurrent->Matrixvalue = xsquareRoot.m;
      mpz_init(xsquareRoot.pcurrent->value);
      mpz_set(xsquareRoot.pcurrent->value, quo1);
      //--------------------
      // output quo1 AND...
      //--------------------
      mpz_neg(quo1, quo1);
      Matrix moutput(zero, one, one, quo1);
      xsquareRoot.m = moutput * xsquareRoot.m;
      mpz_neg(quo1, quo1);
      //-------------------
      // AND... input quo1
      //-------------------
      Matrix minput(quo1, one, one, zero);
      xsquareRoot.m = xsquareRoot.m * minput;
      //------
      // shift
      //------
      mpz_set(*pquo1, x0);
      mpz_set(x0, x1);
      mpz_set(x1, quo1);
      mpz_clear(quo1); mpz_clear(rem);
#ifdef DEBUG
      printf("eval_squareRoot() quo1 = "); mpz_out_str(stdout, 10, quo1); printf("\n");
      printf("eval_squareRoot() return = "); mpz_out_str(stdout, 10, *pquo1); printf("\n");
      xsquareRoot.m.matrixPrint();
#endif
      return pquo1; // that's all folks
   }
}

// calculate  new_x <- int((x + int(m.x))/2) and quo <- int(m.x)
// x, m : input
// quo, new_x, status : output
//---------------------------------------------------------//
void newtonIteration(mpz_t *x, const Matrix & m,
                     mpz_t *quo, mpz_t *new_x, int *status) {
//---------------------------------------------------------//
   mpz_t add, rem; mpz_init(add); mpz_init(rem);

   //     /x\
   // m . | | evaluation
   //     \1/
   Vector v2(*x, one);
   Vector v1 = m * v2;
   v1.VectorSimplify();
#ifdef DEBUG
   printf(" newtonIteration() vector = "); // v1.vectorPrint();
   mpz_out_str(stdout, 10, *v1._geta()); printf(" ");
   mpz_out_str(stdout, 10, *v1._getb()); printf("\n");
#endif
   if (mpz_cmp_ui(*v1._getb(), 0) != 0) {
      mpz_fdiv_q(*quo, *v1._geta(), *v1._getb());
#ifdef DEBUG
      printf("newtonIteration() quo = "); mpz_out_str(stdout, 10, *quo); printf("\n");
#endif
      // new_x <- (x + quo) / 2
      mpz_add(add, *x, *quo);
      mpz_fdiv_qr_ui(*new_x, rem, add, 2);
      *status = 0; // no error
   } else {
      *status = 1; // error
   }
   mpz_clear(add); mpz_clear(rem);
}

// fixed point x = m(x) determination
// m : input
// result, status : output
//-----------------------------------------------------//
void fixedPoint(Matrix & m, mpz_t *result, int *status) {
//-----------------------------------------------------//
   mpz_t old_x, new_x; mpz_init(old_x); mpz_init(new_x);
   mpz_t tmp1, tmp2, tmp3; mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
   mpz_t mx1, mx2; mpz_init(mx1); mpz_init(mx2);
   int count = 0, nb_sol, xstatus;

   // number of solutions determination
   mpz_sub(tmp1, *m._getd(), *m._geta());
   if (mpz_cmp_ui(*m._getc(), 0) == 0) { // if c == 0
      if (mpz_cmp_ui(tmp1, 0) == 0) {
         // d-a == 0
         nb_sol = 0; // 0 solution
      } else {
         // d-a != 0
         nb_sol = 1; // 1 solution
      }
   } else { // case c != 0
      // (d-a)^2 + 4bc evaluation
      mpz_pow_ui(tmp1, tmp1, 2);
      mpz_mul(tmp2, *m._getb(), *m._getc());
      mpz_mul_ui(tmp2, tmp2, 4);
      mpz_add(tmp1, tmp1, tmp2);
      if (mpz_cmp_ui(tmp1, 0) > 0) { // case (d-a)^2 + 4bc  > 0
         nb_sol = 2; // 2 solutions
      } else if (mpz_cmp_ui(tmp1, 0) < 0) { // case (d-a)^2 + 4bc  < 0
         nb_sol = 0; // 0 solution
      } else { // case (d-a)^2 + 4bc  = 0
         nb_sol = 1; // 1 solution
      }
   }
#ifdef DEBUG
   printf("fixedPoint() nb_sol = %d\n", nb_sol);
   printf("fixedPoint() "); m.matrixPrint();
#endif
   mpz_set_ui(tmp1, 1);
rebelotte:
   mpz_set_ui(old_x, 0);
   mpz_set(new_x, tmp1);
rebelotte2:
   while (1) {
#ifdef DEBUG
      printf("fixedPoint() new_x = "); mpz_out_str(stdout, 10, new_x); printf("\n");
#endif
      mpz_set(old_x, new_x); // shift
      // newton iteration : new_x <- (old_x + M(old_x)) / 2 and mx1 <- M(old_x)
      newtonIteration(&old_x, m, &mx1, &new_x, &xstatus);
      if (xstatus == 1) {
         *status = 1; // error propagation
         mpz_clear(old_x); mpz_clear(new_x);
         mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
         mpz_clear(mx1); mpz_clear(mx2);
#ifdef ZDEBUG
         printf("*** error 1 newtonIteration\n");
#endif
         return;
      }
      // mx2 <- M(old_x+1)
      mpz_add_ui(tmp2, old_x, 1);
      newtonIteration(&tmp2, m, &mx2, &tmp3, &xstatus);
      if (xstatus == 1) {
         *status = 1; // error propagation
         mpz_clear(old_x); mpz_clear(new_x);
         mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
         mpz_clear(mx1); mpz_clear(mx2);
#ifdef ZDEBUG
         printf("*** error 2 newtonIteration\n");
#endif
         return;
      }
      // if x = (x + M(x))/2 then x is a fixed point
      if (mpz_cmp(old_x, new_x) == 0) {
         break; // exit from while (1)
      }
#ifdef DEBUG
      printf("fixedPoint() old_x : "); mpz_out_str(stdout, 10, old_x); printf("\n");
      printf("fixedPoint() mx1 : "); mpz_out_str(stdout, 10, mx1); printf("\n");
      printf("fixedPoint() tmp2 : "); mpz_out_str(stdout, 10, tmp2); printf("\n");
      printf("fixedPoint() mx2 : "); mpz_out_str(stdout, 10, mx2); printf("\n");
#endif
      // if x <= M(x) and (x+1) > M(x+1) then x is a fixed point
      if ((mpz_cmp(old_x, mx1) <= 0) && (mpz_cmp(tmp2, mx2) > 0)) {
         mpz_set(new_x, old_x); // fp is old_x
         break; // exit from while (1)
      }
   } // end while (1)
   // negative fixed point elimination
   if (mpz_cmp_ui(new_x, 0) < 0) {
      if (nb_sol == 2) { // find other fixed point
         mpz_mul_ui(tmp1, tmp1, 2);
         count = 0;
#ifdef DEBUG
         printf("fixedPoint() : *2 rebelotte : ");
         mpz_out_str(stdout, 10, tmp1); printf("\n");
#endif
         goto rebelotte;
      } else { // no other fixed point
         *status = 1;
         mpz_clear(old_x); mpz_clear(new_x);
         mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
         mpz_clear(mx1); mpz_clear(mx2);
         return;
      }
   }
   // null fixed point elimination
   if (mpz_cmp_ui(new_x, 0) == 0) {
      count++;
      if (count == 2) {
#ifdef DEBUG
         printf("fixedPoint() fixed point 0 rejected\n");
#endif
         if (nb_sol == 2) { // find other fixed point
           mpz_mul_ui(tmp1, tmp1, 2);
           count = 0;
#ifdef DEBUG
           printf("fixedPoint() : *2 rebelotte : ");
           mpz_out_str(stdout, 10, tmp1); printf("\n");
#endif
           goto rebelotte;
         } else { // no other fixed point
            *status = 1;
            mpz_clear(old_x); mpz_clear(new_x);
            mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
            mpz_clear(mx1); mpz_clear(mx2);
            return;
         }
      } else {
#ifdef DEBUG
         printf("fixedPoint() : rebelotte2\n");
#endif
         goto rebelotte2;
      }
   }
   mpz_set(*result, new_x);
#ifdef DEBUG
   printf("fixedPoint() : "); m.matrixPrint();
   printf("fixedPoint() : "); mpz_out_str(stdout, 10, new_x); printf(" accepted\n");
#endif
   *status = 0;
   mpz_clear(old_x); mpz_clear(new_x);
   mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
   mpz_clear(mx1); mpz_clear(mx2);
}

// squareRootOperation evaluator : square root of normNumber
//--------------------------------------------------//
mpz_t *ContinuedFraction::eval_squareRootOperation() {
//--------------------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::eval_squareRootOperation() %d\n", this);
   xsquareRootOperation.t.tensorPrint();
#endif
   int xoutput;
   Vector v1, v2;
   mpz_t result1, result2; mpz_init(result1); mpz_init(result2);
   mpz_t result3, result4; mpz_init(result3); mpz_init(result4);
   int status;

   // rajouter un test pour voir si x0 vaut 0 !!!
   if (mpz_cmp_ui(xsquareRootOperation.unary->x0, 0) == 0) {
   }

   // if end of squareRootOperation
   //if (mpz_cmp_ui(xsquareRootOperation.unary->x1, 0) == 0) {
   if (mpz_cmp_si(xsquareRootOperation.unary->x1, -77) == 0) {
      //-----------------------------------
      // squareRootOperation -> squareRoot
      //-----------------------------------
#ifdef YDEBUG
      printf("\nsquarerootoperation -> squareroot\n");
#endif
      discriNode = squareRoot;
      Vector v2(xsquareRootOperation.unary->x0, one);
      xsquareRoot.m = xsquareRootOperation.t * v2;
#ifdef DEBUG
      //xsquareRoot.m.matrixPrint();
#endif
      squareRootRat(*xsquareRoot.m._geta(), *xsquareRoot.m._getb(),
                    *xsquareRoot.m._getc(), discriOutput);
      mpz_clear(result1); mpz_clear(result2);
      mpz_clear(result3); mpz_clear(result4);
      return eval(); // recursive call
   // particular case : sqrt([0 x...]) = [0 y...]
   } else if (mpz_cmp_ui(xsquareRootOperation.unary->x0, 0) == 0) {
      //---------
      // y input
      //---------
      Matrix minput2(*(xsquareRootOperation.unary->eval()), // recursive call
                     one, one, zero);
      xsquareRootOperation.t = xsquareRootOperation.t * minput2;
      if ((mpz_cmp_ui(xsquareRootOperation.unary->x0, 1) == 0) &&
          (mpz_cmp_ui(xsquareRootOperation.unary->x1, 0) == 0)) {
         // particular case : sqrt([0 1]) = [1]
         // A MODIFIER !!! sqrt([0 1]) /= [0 1] 
      }
      //---------------------------
      // z output (zero) AND...
      //---------------------------
      Matrix moutput(zero, one, one, zero);
      xsquareRootOperation.t = moutput * xsquareRootOperation.t;
      //--------------------------
      // ...AND x input (zero)
      //--------------------------
      // two tensorTranspose because of x input (and not y input)
      xsquareRootOperation.t.tensorTranspose();
      Matrix minput(zero, one, one, zero);
      xsquareRootOperation.t = xsquareRootOperation.t * minput;
      xsquareRootOperation.t.tensorTranspose();
      mpz_t *pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result
      mpz_set(*pquo1, x0);
      mpz_set(x0, x1);
      mpz_set_ui(x1, 0);
      mpz_clear(result1); mpz_clear(result2);
      mpz_clear(result3); mpz_clear(result4);
      return pquo1; // that's all folks
   } else {
      Vector v3(zero, one), v4(one, zero);
      // y = 0 evaluation
      Matrix m1 = xsquareRootOperation.t * v3;
      // y = oo evaluation
      Matrix m2 = xsquareRootOperation.t * v4;
#ifdef DEBUG
      printf("ContinuedFraction::eval_squareRootOperation() : tensor ");
      xsquareRootOperation.t.tensorPrint();
      printf("ContinuedFraction::eval_squareRootOperation() : (y=0)  righthand ");
      m1.matrixPrint();
      printf("ContinuedFraction::eval_squareRootOperation() : (y=oo) lefthand  ");
      m2.matrixPrint();
#endif
      // fixed point x = M(y=0)(x) determination
      fixedPoint(m1, &result1, &status);
      if (status == 0) {
#ifdef DEBUG
         printf("ContinuedFraction::eval_squareRootOperation() : fixed point 1 : ");
         mpz_out_str(stdout, 10, result1);
         printf("\n");
#endif
         // verification
         // result2 = M(y=oo)(x=fp1=result1) evaluation
         // result3 = M(y=oo)(x=fp1+1) evaluation
         mpz_add_ui(result4, result1, 1);
         Vector v3(result1, one), v4(result4, one);
         v1 = m2 * v3;
         v1.VectorSimplify();
         v2 = m2 * v4;
         v2.VectorSimplify();
#ifdef DEBUG
         printf("v1 = ("); mpz_out_str(stdout, 10, *v1._geta()); printf(" ");
         mpz_out_str(stdout, 10, *v1._getb()); printf(")\n");
         printf("v2 = ("); mpz_out_str(stdout, 10, *v2._geta()); printf(" ");
         mpz_out_str(stdout, 10, *v2._getb()); printf(")\n");
#endif
         if ((mpz_cmp_ui(*v1._getb(), 0) != 0) &&
             (mpz_cmp_ui(*v2._getb(), 0) != 0)) {
            mpz_tdiv_q(result2, *v1._geta(), *v1._getb());
            mpz_tdiv_q(result3, *v2._geta(), *v2._getb());
            // fixed point comparison
            //if (mpz_cmp(result1, result2) == 0) {}
            if ((mpz_cmp(result2, result1) >= 0) &&
                (mpz_cmp(result3, result4) < 0)) {
#ifdef DEBUG
               printf("ContinuedFraction::eval_squareRootOperation() : output/input\n");
#endif
               xoutput = 1; // output
            } else {
#ifdef DEBUG
               printf("ContinuedFraction::eval_squareRootOperation() : fp1 != fp2 : y input\n");
#endif
               xoutput = 0; // input
            }
         } else {
#ifdef DEBUG
            printf("ContinuedFraction::eval_squareRootOperation() : v1.b == 0 or v2.b == 0\n");
#endif
            xoutput = 0; // input
         }
      } else {
#ifdef ZDEBUG
         printf("ContinuedFraction::eval_squareRootOperation() : no fixed point : y input\n");
#endif
         xoutput = 0; // input
      }
   }
   if (xoutput) {
      if (discriOutput == decimalOutput) {
         printf("impossible error...\n");
         exit(-1);
      }
      //---------------------------
      // z output (result1) AND...
      //---------------------------
      mpz_neg(result1, result1);
      Matrix moutput(zero, one, one, result1);
      xsquareRootOperation.t = moutput * xsquareRootOperation.t;
#ifdef TDEBUG
      printf("*** output "); mpz_out_str(stdout, 10, result1); printf("\n");
      xsquareRootOperation.t.tensorPrint();
#endif
      mpz_neg(result1, result1);
      //--------------------------
      // ...AND x input (result1)
      //--------------------------
      // two tensorTranspose because of x input (and not y input)
      xsquareRootOperation.t.tensorTranspose();
      Matrix minput(result1, one, one, zero);
      xsquareRootOperation.t = xsquareRootOperation.t * minput;
      xsquareRootOperation.t.tensorTranspose();
#ifdef TDEBUG
      printf("*** input "); mpz_out_str(stdout, 10, result1); printf("\n");
      xsquareRootOperation.t.tensorPrint();
#endif
      mpz_t *pquo1 = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*pquo1); // for the result
      mpz_set(*pquo1, x0);
      mpz_set(x0, x1);
      mpz_set(x1, result1);
      mpz_clear(result1); mpz_clear(result2);
      mpz_clear(result3); mpz_clear(result4);
#ifdef DEBUG
      printf("ContinuedFraction::eval_squareRootOperation() : fin ");
      mpz_out_str(stdout, 10, *pquo1);
      // printf(" "); mpz_out_str(stdout, 10, result1);
      printf("\n");
#endif
      return pquo1; // that's all folks
   } else {
      //---------
      // y input
      //---------
      Matrix minput(*(xsquareRootOperation.unary->eval()),
                    one, one, zero); // recursive call
      // free(result of unary->eval); ???
      xsquareRootOperation.t = xsquareRootOperation.t * minput;
      mpz_clear(result1); mpz_clear(result2);
      mpz_clear(result3); mpz_clear(result4);
      return eval(); // recursive call
   }
}

//--------------------------------------//
Vector ContinuedFraction::eval_glouton() {
//--------------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::eval_glouton() %d\n", this);
   xunaryOperation.m.matrixPrint();
#endif

   // end of unaryOperation
   //if (mpz_cmp_ui(xunaryOperation.unary->x1, 0) == 0) {
   if (mpz_cmp_si(xunaryOperation.unary->x1, -77) == 0) {
      //----------------------------
      // unaryOperation -> rational
      //----------------------------
#ifdef YDEBUG
      printf("\nunaryoperation -> rational\n");
#endif
#ifdef DEBUG
      printf("unaryOperation -> rational ");
      mpz_out_str(stdout, 10, xunaryOperation.unary->x0);
      printf("\n");
      //printf("("); mpz_out_str(stdout, 10, *v._geta()); printf(" ");
      //             mpz_out_str(stdout, 10, *v._getc()); printf(")\n");
#endif
      discriNode = rational; // end of eval_glouton()
      Vector v2(xunaryOperation.unary->x0, one);
      v = xunaryOperation.m * v2;
      v.VectorSimplify();
      return v;
   } else {
      //---------
      // input a
      //---------
      Matrix minput(*(xunaryOperation.unary->eval()), one, one, zero); // recursive call
      // free(result of unary->eval); ???
      xunaryOperation.m = xunaryOperation.m * minput;
#ifdef DEBUG
      printf("eval_glouton : input ");
      //matrixPrint(xunaryOperation.m);
#endif
      Vector vresult(*xunaryOperation.m._geta(), *xunaryOperation.m._getc());
      return vresult;
   }
}

//---------------------
// remarquable numbers
//---------------------

// coth(1) = (exp(1)+exp(-1))/(exp(1)-exp(-1)) = [1 3 5 7 9 11...]
//-----------------------//
mpz_t *coth1(mpz_t state) { // for normNumber
//-----------------------//
   mpz_t *result = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*result); // for the result

   mpz_mul_ui(*result, state, 2);
   mpz_add_ui(*result, *result, 1); // 2*state+1
   return result;
}

// e = neperien logarithms base = [2 1 2 1 1 4 1 1 6 1 1 8...]
//-------------------//
mpz_t *e(mpz_t state) { // for normNumber
//-------------------//
   mpz_t tmp; mpz_init(tmp);
   mpz_t *result = (mpz_t *)malloc(sizeof(mpz_t)); mpz_init(*result); // for the result

   mpz_add_ui(tmp, state, 1);
   mpz_mul_ui(tmp, tmp, 2);
   if (mpz_cmp_ui(state, 0) == 0) { // state == 0
      mpz_set_ui(*result, 2);
   } else if (mpz_divisible_ui_p(tmp, 3)) { // (state+1)%3 == 0
      mpz_divexact_ui(*result, tmp, 3); // 2*(state+1)/3
   } else {
      mpz_set_ui(*result, 1);
   }
   mpz_clear(tmp);
   return result;
}

//              4
//  pi  =  0 + -----------------------------
//                   1
//              1 + ------------------------
//                        4
//                   3 + -------------------
//                              9
//                        5 + --------------
//                                  16
//                             7 + ---------
//                                  9 + ...
//--------------------//
Vector pi(mpz_t state) { // for unnormNumber
//--------------------//
   Vector result;
   mpz_t tmp; mpz_init(tmp);

   if (mpz_cmp_ui(state, 0) == 0) { // state == 0
      mpz_set_ui(*result._geta(), 0);         // a = 0
      mpz_set_ui(*result._getb(), 4);         // b = 4
   } else {
      mpz_mul_ui(tmp, state, 2);
      mpz_sub_ui(*result._geta(), tmp, 1);    // a = 2*state-1
      mpz_mul(*result._getb(), state, state); // b = state*state
   }
   mpz_clear(tmp);
   return result;
}

//      n    /0  n\   /0    n\   /0    n\   /0     n\   /0     n\
// atan(-) = |    | . |      | . |      | . |       | . |       | ....
//      d    \n  d/   \4n  3d/   \9n  5d/   \16n  7d/   \25n  9d/
//            k=0       k=1        k=2        k=3         k=4
//-------------------------------------------//
Matrix atanRat(mpz_t state, mpz_t n, mpz_t d) { // for generalizedUnnormNumber
//-------------------------------------------//
   Matrix mresult;
   mpz_t tmp; mpz_init(tmp);

   mpz_set_ui(*mresult._geta(), 0);   // a = 0
   mpz_set(*mresult._getb(), n);      // b = n
   mpz_add_ui(tmp, state, 1);
   mpz_mul(tmp, tmp, tmp);
   mpz_mul(*mresult._getc(), n, tmp); // c = n*(state+1)*(state+1)
   mpz_mul_ui(tmp, state, 2);
   mpz_add_ui(tmp, tmp, 1);
   mpz_mul(*mresult._getd(), d, tmp); // d = d*(2*state+1)
   mpz_clear(tmp);
   return mresult;
}

//     n    /0  n\   /0   n\   /0   n\   /0   n\   /0   n\   /0    n\
// tan(-) = |    | . |     | . |     | . |     | . |     | . |      | ...
//     d    \n  d/   \n -3d/   \n  5d/   \n -7d/   \n  9d/   \n -11d/
//           k=0       k=1       k=2       k=3       k=4       k=5
//-------------------------------------------//
Matrix tanRat(mpz_t state, mpz_t n, mpz_t d) { // for generalizedUnnormNumber
//-------------------------------------------//
   Matrix mresult;
   mpz_t tmp; mpz_init(tmp);

   mpz_set_ui(*mresult._geta(), 0);   // a = 0
   mpz_set(*mresult._getb(), n);      // b = n
   mpz_set(*mresult._getc(), n);      // c = n

   mpz_mul_ui(tmp, state, 2);
   mpz_add_ui(tmp, tmp, 1);
   mpz_mul(*mresult._getd(), d, tmp); // d = d*(2*state+1)
   if (mpz_divisible_ui_p(state, 2)) { // state%2 == 0 (state even)
   } else { // state odd
      mpz_neg(*mresult._getd(), *mresult._getd()); // d = - d*(2*state+1)
   }
   mpz_clear(tmp);
   return mresult;
}

//      n          p    /0  p\   /0  1\   /0   p\   /0  2\   /0   p\   /0  3\
// log (-) = log(1+-) = |    | . |    | . |     | . |    | . |     | . |    | ...
//      d          q    \p  q/   \1  2/   \p  3q/   \2  2/   \p  5q/   \3  2/
//                       k=0      k=1       k=2      k=3       k=4      k=5
//------------------------------------------//
Matrix logRat(mpz_t state, mpz_t n, mpz_t d) { // for generalizedUnnormNumber
//------------------------------------------//
   Matrix mresult;
   mpz_t tmp; mpz_init(tmp);
   mpz_t p, q; mpz_init(p); mpz_init(q);

   mpz_sub(p, n, d); // p = n - d
   mpz_set(q, d);    // q = d
   if (mpz_divisible_ui_p(state, 2)) { // state%2 == 0 (state even)
      mpz_set_ui(*mresult._geta(), 0);           // a = 0
      mpz_set(*mresult._getb(), p);              // b = p
      mpz_set(*mresult._getc(), p);              // c = p
      mpz_add_ui(tmp, state, 1);
      mpz_mul(*mresult._getd(), q, tmp);         // d = (state+1)*q
   } else { // state odd
      mpz_set_ui(*mresult._geta(), 0);           // a = 0
      mpz_add_ui(tmp, state, 1);
      mpz_divexact_ui(*mresult._getb(), tmp, 2); // b = (state+1)/2
      mpz_divexact_ui(*mresult._getc(), tmp, 2); // c = (state+1)/2
      mpz_set_ui(*mresult._getd(), 2);           // d = 2
   }
   mpz_clear(tmp);
   mpz_clear(p); mpz_clear(q);
   return mresult;
}

//     n    /2d+n  n\   /6d  n\   /10d  n\   /14d  n\
// exp(-) = |       | . |     | . |      | . |      | ...
//     d    \2d-n  n/   \n   0/   \n    0/   \n    0/
//             k=0        k=1        k=2        k=3
//------------------------------------------//
Matrix expRat(mpz_t state, mpz_t n, mpz_t d) { // for generalizedUnnormNumber
//------------------------------------------//
   Matrix mresult;
   mpz_t tmp; mpz_init(tmp);

   if (mpz_cmp_ui(state, 0) == 0) { // state == 0
      mpz_mul_ui(tmp, d, 2);
      mpz_add(*mresult._geta(), tmp, n);       // a = 2*d+n
      mpz_set(*mresult._getb(), n);            // b = n
      mpz_sub(*mresult._getc(), tmp, n);       // c = 2*d-n
      mpz_set(*mresult._getd(), n);            // d = n
   } else {
      mpz_mul_ui(tmp, state, 4);
      mpz_add_ui(tmp, tmp, 2);
      mpz_mul(*mresult._geta(), d, tmp);       // a = (4*state+2)*d
      mpz_set(*mresult._getb(), n);            // b = n
      mpz_set(*mresult._getc(), n);            // c = n
      mpz_set_ui(*mresult._getd(), 0);         // d = 0
   }
   mpz_clear(tmp);
   return mresult;
}

//      case 14: // a normalized number node coth(1)
//               //----------------------------------
//         cf = ContinuedFraction(&coth1, discriOutput);
//      case 16: // squareRootOperation : coth(1/2)
//               //---------------------------------
//         cf1 = ContinuedFraction(&coth1, ContinuedFraction::cfOutput);
//         cf1.x0Andx1Init();
//         cf = cf.cothHalfOne(cf1, discriOutput);
//      case 17: // squareRootOperation : sqrt(coth(1))
//               //-------------------------------------
//         cf1 = ContinuedFraction(&coth1, ContinuedFraction::cfOutput);
//         cf1.x0Andx1Init();
//         cf = cf.squareRootOp(cf1, discriOutput);

//--------------------------------//
void ContinuedFraction::evalMain() {
//--------------------------------//
#ifdef DEBUG
   printf("ContinuedFraction::evalMain()\n");
#endif
   //x0Andx1Init();
   int i = 0; // current length output
   mpz_t xresult; mpz_init(xresult);

   printf("[");

   // FIRST output
   //--------------
   if ((discriOutput == ContinuedFraction::cfOutput) ||
       (discriOutput == ContinuedFraction::decimalOutput)) {
      //-----------------------------------------
      // first call of ContinuedFraction::eval()
      //-----------------------------------------
      mpz_set(xresult, *(eval()));
      //-----------------------------------------
      // free(result of eval()); ???
      // ATTENTION !!! the FIRST result must be null...
      mpz_out_str(stdout, 10, xresult);
   } else if (discriOutput == ContinuedFraction::rationalOutput) {
      //-------------------------------------------------
      // first call of ContinuedFraction::eval_glouton()
      //-------------------------------------------------
      Vector xvector = eval_glouton();
      //-------------------------------------------------
#if 0
      printf(" eval_glouton() vector = ");
      mpz_out_str(stdout, 10, *xvector._geta()); printf(" ");
      mpz_out_str(stdout, 10, *xvector._getb()); printf("\n");
#endif
   }
   i++;

   if (mpz_cmp_ui(scale, 0) > 0) { // scale : length of cf output
   // SECOND, THIRD... outputs
   //--------------------------
   // case cfOutput
   //---------------
   if (discriOutput == ContinuedFraction::cfOutput) {
      while (1) {
         //-----------------------------------------
         // next calls of ContinuedFraction::eval()
         //-----------------------------------------
         mpz_set(xresult, *(eval()));
         //-----------------------------------------
         // free(result of eval()); ???
         //if (mpz_cmp_ui(xresult, 0) == 0) {
         if (mpz_cmp_si(xresult, -77) == 0) {
            break; // exit while (1)
         } else if (mpz_cmp_ui(scale, i) == 0) { // scale : length of cf output
            printf(" ");
            mpz_out_str(stdout, 10, xresult);
            printf("...");
            break; // exit while (1)
         } else {
            printf(" ");
            mpz_out_str(stdout, 10, xresult);
            i++;
         }
      }
   // case decimalOutput
   //--------------------
   } else if (discriOutput == ContinuedFraction::decimalOutput) {
      if (mpz_cmp_ui(xresult, 0) < 0) {
         printf("+0");
      }
      printf(".");
      while (1) {
         //-----------------------------------------
         // next calls of ContinuedFraction::eval()
         //-----------------------------------------
         mpz_set(xresult, *(eval()));
         //-----------------------------------------
         // free(result of eval()); ???
         // ATTENTION !!! the results must be null in case of decimal output...
         if (mpz_cmp_ui(scale, i) == 0) { // scale : length of decimal output
            mpz_out_str(stdout, 10, xresult);
            printf("...");
            break; // exit while (1)
         } else {
            mpz_out_str(stdout, 10, xresult);
            i++;
         }
      }
   } else if (discriOutput == ContinuedFraction::rationalOutput) {
      while (1) {
         Vector xvector = eval_glouton();
         i++;
         if (mpz_cmp_ui(scale, i) == 0) { // scale : length of cf output
            mpz_out_str(stdout, 10, *xvector._geta());
            printf("/");
            mpz_out_str(stdout, 10, *xvector._getb());
            break; // exit while (1)
         }
      }
   }
   }

   printf("]\n");
   mpz_clear(xresult);
}

extern int yyparse();
extern void scan_string(char *); // cf cf_arith2.l
extern void delete_string(); // cf cf_arith2.l

jmp_buf env; // setjmp/longjmp

//--------------------------//
void signal_handler(int sig) {
//--------------------------//
   switch (sig) {
      case SIGINT:
         longjmp(env, sig);
         // never reached
      case SIGALRM:
         longjmp(env, sig);
         // never reached
      default:
         exit(sig);
   }
}

//-----------------------------//
int main(int argc, char **argv) {
//-----------------------------//
   FILE *in, *out;
   char *infile;
   char *line;

   int returned_from_longjmp;
   //void *signal_alarm_save; // ancienne déclaration provoquant une erreur
   void (*signal_alarm_save) (int);

   if (argc > 1) {
      infile = argv[1];
      in = freopen(infile, "r", stdin);
      if (in == NULL) {
         fprintf(stderr, "cannot open %s\n", infile);
         exit(ENOENT);
      }
   } else {
      in = stdin;
   }
   mpz_init(zero); mpz_set_ui(zero, 0);
   mpz_init(one);  mpz_set_ui(one,  1);
   mpz_init(ten);  mpz_set_ui(ten,  10);
   mpz_init(scale); mpz_set_ui(scale, 50); // default output length
   discriOutput = ContinuedFraction::cfOutput; // default output type

   // setjmp processing
   if ((returned_from_longjmp = setjmp(env)) != 0) {
      // ici, il faudrait faire tous les free de evalMain() : tous les ~
      switch (returned_from_longjmp) {
         case SIGINT:
            printf("*** longjumped from interrupt %d\n", SIGINT);
            break;
         case SIGALRM:
            printf("*** longjumped from alarm %d\n", SIGALRM);
            break;
      }
   }

   while (1) {
      line = readline(">"); // output prompt + readline/history library
      if (line == NULL) {
          break;
      }
      if (line && *line) {
         add_history(line); // readline/history library
      }
      //if (*line != NULL) { // a cause de NULL used in arithmetic
      if (*line) {
         // signal processing
         //(void)signal(SIGINT, signal_handler);
         signal_alarm_save = signal(SIGALRM, signal_handler);
         alarm(30); // in seconds

         scan_string(line); // lex scan of the read line (cf cf_arith2.l)
         yyparse();
         delete_string();

         alarm(0); // reset alarm
      }
      free(line);
   }
   fclose(in);
}

extern int yyerror(char *);

struct symbolTable symbolTable[NBSYM];

// find symbol in the table if it exists...
// and adds it if it doesn't...
struct symbolTable *findSymbol(char *s) {
   for (struct symbolTable *sp = symbolTable; sp < &symbolTable[NBSYM]; sp++) {
      if (sp->name && !strcmp(sp->name, s)) // already present ?
         return sp;
      if (!sp->name) { // free entry ?
         // new entry
         sp->name = strdup(s); // implicit malloc
         // initialization : value <- ZERO
         sp->value = new ContinuedFraction(0, 1, ContinuedFraction::cfOutput);
         return sp;
      }
   }
   yyerror("symbol table full...\n");
   // exit(1);
   return NULL;
}

/* (akozar) MOVED TO cf_arith.l:
int yywrap(void) { return 1;}
*/
