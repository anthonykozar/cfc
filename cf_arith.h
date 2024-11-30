//--------------------------------------------------------------------//
// continued fractions arithmetic with Gosper's algorithm (cf HAKMEM) //
//--------------------------------------------------------------------//
// the multiple precision GNU library GMP is required (cf http://www.swox.com/gmp/)

#include "gmp.h" // APIs prototypes of the bignum library

//#define DEBUG

//-------------------------//
// constants : 0, 1 and 10 //
//-------------------------//
extern mpz_t zero, one, ten;
extern mpz_t scale; // length of cf output

class Vector; // forward declaration
class Matrix; // forward declaration
class Tensor; // forward declaration
class ContinuedFraction; // forward declaration

//----------//
class Vector {
//----------//
   private:
      //---------------------------------------//
      //                            a     /a\  //
      //  vector structure :  v  =  -  =  | |  //
      //                            b     \b/  //
      //---------------------------------------//
      mpz_t a, b;
   public:
      //------------------------------
      // constructors and destructors
      //------------------------------
      Vector();
      Vector(mpz_t, mpz_t);
      Vector(int, int);
      ~Vector();
      mpz_t * _geta() { return &a; }
      mpz_t * _getb() { return &b; }
      Vector(const Vector &); // initialization
      Vector & operator=(const Vector &); // assignment
      void VectorSimplify();
      //--------------
      // v  =  m1 * v
      //--------------
      friend const Vector operator*(const Matrix & m1, const Vector & v); // forward declaration
};

//----------//
class Matrix {
//----------//
   private:
      //----------------------------------------------//
      //                            ax + b     /a b\  //
      //  matrix structure :  m  =  ------  =  |   |  //
      //                            cx + d     \c d/  //
      //----------------------------------------------//
      mpz_t a, b, c, d;
      void init(); // private initialization
   public:
      //------------------------------
      // constructors and destructors
      //------------------------------
      Matrix();
      Matrix(mpz_t, mpz_t, mpz_t, mpz_t);
      Matrix(int, int, int, int);
      Matrix(Vector & v1, Vector  &v2);
      ~Matrix();
      mpz_t * _geta() { return &a; }
      mpz_t * _getb() { return &b; }
      mpz_t * _getc() { return &c; } // attention getc() collisionne
      mpz_t * _getd() { return &d; }
      Matrix(const Matrix &); // initialization
      const Matrix & operator=(const Matrix &); // assignment
      //-----------
      // m  +=  m1
      //-----------
      const Matrix & operator+=(const Matrix &);
      //---------------
      // m  =  m1 + m2
      //---------------
      friend const Matrix operator+(const Matrix &, const Matrix &);
      //-----------
      // m  *=  m1
      //-----------
      const Matrix & operator*=(const Matrix &);
      //---------------
      // m  =  m1 * m2
      //---------------
      friend const Matrix operator*(const Matrix &, const Matrix &);
      //--------------
      // m  =  i * m1
      //--------------
      friend const Matrix operator*(const mpz_t, const Matrix &);
      void matrixPrint(); // print (a b c d)
      void matrixTranspose();
      //----------------------------------------------------
      //                 /ma mb\   /va\     /ma*va + mb*vb\
      // v  =  m * v  =  |     | * |  |  =  |             |
      //                 \mc md/   \vb/     \mc*va + md*vb/
      //----------------------------------------------------
      friend const Vector operator*(const Matrix & m1, const Vector & v);
      friend const Tensor operator*(const Matrix &, const Tensor &); // forward declaration
      friend const Tensor operator*(const Tensor & t, const Matrix & m); // forward declaration
};

//----------//
class Tensor {
//----------//
   private:
      //---------------------------------------------------------------//
      //  tensor structure :                                           //
      //                                                               //
      //        axy + bx + cy + d     / /a b\  /e f\ \     /        \  //
      //  t  =  -----------------  =  | |   |, |   | |  =  | m1, m2 |  //
      //        exy + fx + gy + h     \ \c d/  \g h/ /     \        /  //
      //                                                               //
      //---------------------------------------------------------------//
      Matrix matrix_1, matrix_2;
   public:
      //------------------------------
      // constructors and destructors
      //------------------------------
      Tensor();
      Tensor(mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t);
      Tensor(int, int, int, int, int, int, int, int);
      Tensor(Matrix &, Matrix &);
      ~Tensor();
      mpz_t * _geta() { return matrix_1._geta(); }
      mpz_t * _getb() { return matrix_1._getb(); }
      mpz_t * _getc() { return matrix_1._getc(); }
      mpz_t * _getd() { return matrix_1._getd(); }
      mpz_t * _gete() { return matrix_2._geta(); }
      mpz_t * _getf() { return matrix_2._getb(); }
      mpz_t * _getg() { return matrix_2._getc(); }
      mpz_t * _geth() { return matrix_2._getd(); }
      Matrix _getm1() const { return matrix_1; }
      Matrix _getm2() const { return matrix_2; }
      Tensor(const Tensor &); // initialization
      Tensor & operator=(const Tensor &); // assignment
      void tensorPrint(); // print (a b c d) / (e f g h)
      void tensorTranspose();
      //-----------------------------------------------------
      //                 //      \   /v1\\T     /         \T
      // m  =  t * v  =  ||m1, m2| * |  ||   =  |m1.v m2.v|
      //                 \\      /   \v2//      \         /
      //-----------------------------------------------------
      friend const Matrix operator*(const Tensor &, const Vector &);
      //-------------------------------------------------------------
      //                 /a b\   /      \     /                    \
      // t  =  m * t  =  |   | * |m1, m2|  =  |a*m1+b*m2, c*m1+d*m2|
      //                 \c d/   \      /     \                    /
      //-------------------------------------------------------------
      friend const Tensor operator*(const Matrix &, const Tensor &);
      //-----------------------------------------------
      //                 /      \         /          \
      // t  =  t * m  =  |m1, m2| * m  =  |m1.m, m2.m|
      //                 \      /         \          /
      //-----------------------------------------------
      friend const Tensor operator*(const Tensor & t, const Matrix & m);
};

// functions returning next output
//---------------------------------
typedef mpz_t * (*normNumberFunction)             (mpz_t state);
typedef Vector  (*unnormNumberFunction)           (mpz_t state);
typedef Matrix  (*generalizedUnnormNumberFunction)(mpz_t state, mpz_t n, mpz_t d);

// linked list for finite or periodic continued fraction
struct chain {
   mpz_t value;
   struct chain *next;
};

struct rationalcf {
   mpz_t num, den;
};

// linked list for periodicity detection of rational square root
struct chain2 {
   Matrix Matrixvalue;
   mpz_t value;
   struct chain2 *next;
};

//---------------------//
class ContinuedFraction {
//---------------------//
   public:
      // type of node's tree
      //---------------------
      enum enumNode {
         rational,                // a vector to normalize
         normNumber,              // digit with function returning a number
         normNumber2,             // finite or periodical normNumber (linked list)
         unnormNumber,            // digit with function returning a vector to normalize
         generalizedUnnormNumber, // digit with matrix to normalize
         operation,               // tensor with two pointers
         unaryOperation,          // matrix with one pointer
         squareRoot,
         squareRootOperation
      };
      // type of output
      //----------------
      enum enumOutput {
         cfOutput,
         decimalOutput,
         rationalOutput           // accept all input, return a rational (Vector)
      };
      //------------------------------
      // constructors and destructors
      //------------------------------
      ContinuedFraction();
      ~ContinuedFraction();
      // rational
      ContinuedFraction(int, int, enumOutput);
      ContinuedFraction(mpz_t, mpz_t, enumOutput);
      // normNumber
      ContinuedFraction(normNumberFunction, enumOutput);
      // normNumber2
      ContinuedFraction(struct chain *, enumOutput);
      // unnormNumber
      ContinuedFraction(unnormNumberFunction, enumOutput);
      // generalizedUnnormNumber
      ContinuedFraction(generalizedUnnormNumberFunction, mpz_t, mpz_t, enumOutput);
      // unaryOperation
      ContinuedFraction(Matrix &, ContinuedFraction &, enumOutput);
      ContinuedFraction(Matrix & m, ContinuedFraction &);
      // operation
      ContinuedFraction(Tensor &, ContinuedFraction &, ContinuedFraction &, enumOutput);
      ContinuedFraction plus(ContinuedFraction &, ContinuedFraction &, enumOutput);
      ContinuedFraction substract(ContinuedFraction &, ContinuedFraction &, enumOutput);
      ContinuedFraction multiply(ContinuedFraction &, ContinuedFraction &, enumOutput);
      ContinuedFraction divide(ContinuedFraction &, ContinuedFraction &, enumOutput);
      ContinuedFraction squareRootRat(mpz_t, mpz_t, enumOutput);
      ContinuedFraction squareRootRat(mpz_t, mpz_t, mpz_t, enumOutput);
      ContinuedFraction squareRootOp(ContinuedFraction &, enumOutput);
      ContinuedFraction cothHalfOne(ContinuedFraction &, enumOutput);
      void x0Andx1Init();
      // generic evaluator
      //-------------------
      void evalMain();
      // rational evaluator
      //--------------------
      Vector eval_glouton();
   private:
      // private variables
      //-------------------
      enumNode   discriNode;   // rational .. squareRootOperation
      enumOutput discriOutput; // cfOutput decimalOutput
      mpz_t x0; // next output
      mpz_t x1; // next next output
      struct chain *partialEval; // [a ... e] + rest
      // private body struct of nodes
      //------------------------------
      Vector v; // for rational
      struct {
         normNumberFunction function;
         mpz_t state; // function state
      } xnormNumber;
      struct {
         struct chain *pchain;
      } xnormNumber2;
      struct {
         Matrix m;
         unnormNumberFunction function;
         mpz_t state; // function state
      } xunnormNumber;
      struct {
         Matrix m;
         mpz_t n, d;
         generalizedUnnormNumberFunction function;
         mpz_t state; // function state
      } xgeneralizedUnnormNumber;
      struct {
         Tensor t;
         ContinuedFraction *left, *right;
      } xoperation;
      struct {
         Matrix m;
         ContinuedFraction *unary;
      } xunaryOperation;
      struct {
         Matrix m;
         mpz_t sqrt_det;
         struct chain2 *pbegin; // for periodicity detection
         struct chain2 *pcurrent; // for periodicity detection
      } xsquareRoot;
      struct {
         Tensor t;
         ContinuedFraction *unary;
      } xsquareRootOperation;
      // private evaluators
      //--------------------
      mpz_t *eval_rational();
      mpz_t *eval_normNumber();
      mpz_t *eval_normNumber2();
      mpz_t *eval_unnormNumber();
      mpz_t *eval_generalizedUnnormNumber();
      mpz_t *eval_unaryOperation();
      mpz_t *eval_operation();
      mpz_t *eval_squareRoot();
      mpz_t *eval_squareRootOperation();
      mpz_t *eval();
};

extern void evalMain(ContinuedFraction *);
extern Vector pi(mpz_t);
extern mpz_t *e(mpz_t);
extern Matrix atanRat(mpz_t, mpz_t, mpz_t);
extern Matrix tanRat(mpz_t, mpz_t, mpz_t);
extern Matrix logRat(mpz_t, mpz_t, mpz_t);
extern Matrix expRat(mpz_t, mpz_t, mpz_t);
extern mpz_t scale; // length of cf output
extern ContinuedFraction::enumOutput discriOutput; // decimal/continued fraction output

#define NBSYM 20 // max number of symbols

struct symbolTable {
   char *name;
   ContinuedFraction *value;
};

struct symbolTable *findSymbol(char *);

extern jmp_buf env; // setjmp/longjmp

