// Program copyright Carey Bloodworth July 16, 2001  Carey@Bloodworth.org
// Uses the GNU GMP math package.
// Maybe be used, modified, and extended provided proper credit is given
// to me for the program and to Robert William ('Bill') Gosper for the
// continued fraction methods and for the "Matrix Tower" idea.

// Continued fractions based on Robert Gosper's "Matrix Tower" idea.
// gcc -Wall gosper.c -o gosper.exe -lgmp

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "gmp.h"

// Do standard binary splitting or Gosper's tower method?
//#define USE_BINARY_SPLITTING 1
#define USE_MATRIX_TOWER 1

unsigned _stklen = 1048576*16;
// Set the compiler's stack length.
// Partially for the recursion but mostly for GMP, just in case.

#define NAT_LOG_PI 1  // pi calculation
//#define NAT_LOG 1     // e calculation
//#define SQRT2 1       // sqrt(2)
//#define ARCCOTAN 239  // atn(1/239)
//#define ARCCOTAN 5    // atn(1/5)

typedef struct {
   mpz_t NewNum,   OldNum;
   mpz_t NewDenom, OldDenom;
} Point;

void InitPoint(Point *P) {
   mpz_init(P->NewNum);
   mpz_init(P->NewDenom);
   mpz_init(P->OldNum);
   mpz_init(P->OldDenom);
}

void DeInitPoint(Point *P) {
   mpz_clear(P->NewNum);
   mpz_clear(P->NewDenom);
   mpz_clear(P->OldNum);
   mpz_clear(P->OldDenom);
}

void PrintPoint(char *Str, Point *P) {
   if (Str)
      printf("%s", Str);
   mpz_out_str(stdout, 10, P->NewNum); printf("\n");
   mpz_out_str(stdout, 10, P->NewDenom); printf("\n");
   mpz_out_str(stdout, 10, P->OldNum); printf("\n");
   mpz_out_str(stdout, 10, P->OldDenom); printf("\n");
}

void BuildPoint(Point *Tip, int TermNo) {
   // Terms start with '0', the integer part.  By default, it's 1:0
   // which means 0 integer part.
   mpz_set_ui(Tip->OldNum, 0);
   mpz_set_ui(Tip->OldDenom, 1);
   // For term zero, the integer part.  Default of zero.
   if (TermNo == 0) {
      mpz_set_ui(Tip->NewNum, 1);
      mpz_set_ui(Tip->NewDenom, 0);  // Set this if leading integer exists.
   }
#if defined(NAT_LOG_PI)
   // 4/pi=1+1^2:3+2^2:5+3^2:7+...
   if (TermNo == 0) {
      mpz_set_ui(Tip->NewDenom, 1);
   } else {
      mpz_set_ui(Tip->NewNum, TermNo);
      mpz_mul(Tip->NewNum, Tip->NewNum, Tip->NewNum);
      mpz_set_ui(Tip->NewDenom, TermNo);
      mpz_mul_ui(Tip->NewDenom, Tip->NewDenom, 2);
      mpz_add_ui(Tip->NewDenom, Tip->NewDenom, 1);
   }
#elif defined(NAT_LOG)
   // e=2+1:1+1:2+1:1+1:1+1:4+1:1+1:1+1:6+....
   if (TermNo == 0) {
      mpz_set_ui(Tip->NewDenom, 2);
   } else {
      mpz_set_ui(Tip->NewNum, 1);
      mpz_set_ui(Tip->NewDenom, 1);
      if ((TermNo % 3) == 2) {
         mpz_set_ui(Tip->NewDenom, ((TermNo/3)*2)+2);
      }
   }
#elif defined(SQRT2)
   // sqrt(2)=1+1:2+1:2+1:2+...
   if (TermNo == 0) {
      mpz_set_ui(Tip->NewDenom, 1);
   } else {
      mpz_set_ui(Tip->NewNum, 1);
      mpz_set_ui(Tip->NewDenom, 2);
   }
#elif defined(ARCCOTAN)
/* Uses an unusual arctan

 ARCTAN(X)                   1
 ---------  = ---------------------------------
     X                   (1*X)**2
              1 + -----------------------------
                            (2*X)**2
                  3 + -------------------------
                               (3*X)**2
                      5 + ---------------------
                                  .
                                   .
                                    .
                                     (N*X)**2
                         (2*N-1) + ------------
                                     (2*N+1)
                                        .
                                         .
                                          .
*/
// 0+1:1+1^2:3+2^2:5+3^2:7+....
/*
** This is a little hard to do because the 'X' is really 1/y
** because we actually want arccotan rather than arctan, and
** we are trying to build it up from the beginning.  If we were
** doing it from the end, it'd be easy.  The trick is to realize
** that the X (ie: 1/y) will actually end up just multiplying the
** denominator by 'y' *and* multiply the remaining numerators also
** by 'y'.  By doing the numerator, we are actually doing it from
** the *previous* term.  It's just a little late because at that
** point we didn't have this term to do it yet.
*/
#define ARCCOTAN_SQR (ARCCOTAN*ARCCOTAN)
   if (TermNo == 0) {
   } else if (TermNo == 1) {
      mpz_set_ui(Tip->NewNum, 1);
      mpz_set_ui(Tip->NewDenom, 1);
   } else {
      mpz_set_ui(Tip->NewNum, TermNo-1);
      mpz_mul(Tip->NewNum, Tip->NewNum, Tip->NewNum);
      mpz_set_ui(Tip->NewDenom, (TermNo*2)-1);
      mpz_mul_ui(Tip->NewDenom, Tip->NewDenom, ARCCOTAN_SQR);
      if (TermNo != 2)
         mpz_mul_ui(Tip->NewNum, Tip->NewNum, ARCCOTAN_SQR);
   }
#else
#error Undefined type of continued fraction to compute.
#endif
}

void MergePoints(Point *Base, Point *Tip) {
   // Base=Base + Tip
  mpz_t Temp1, Temp2, Temp3, Temp4;
  mpz_t Temp5, Temp6, Temp7, Temp8;

  mpz_init(Temp1);mpz_init(Temp2);mpz_init(Temp3);mpz_init(Temp4);
  mpz_init(Temp5);mpz_init(Temp6);mpz_init(Temp7);mpz_init(Temp8);

  mpz_mul(Temp1, Base->NewNum,   Tip->NewDenom);
  mpz_mul(Temp2, Base->OldNum,   Tip->NewNum);
  mpz_mul(Temp3, Base->NewDenom, Tip->NewDenom);
  mpz_mul(Temp4, Base->OldDenom, Tip->NewNum);

  mpz_mul(Temp5, Base->NewNum,   Tip->OldDenom);
  mpz_mul(Temp6, Base->OldNum,   Tip->OldNum);
  mpz_mul(Temp7, Base->NewDenom, Tip->OldDenom);
  mpz_mul(Temp8, Base->OldDenom, Tip->OldNum);

  mpz_add(Base->NewNum  , Temp1, Temp2);
  mpz_add(Base->NewDenom, Temp3, Temp4);
  mpz_add(Base->OldNum  , Temp5, Temp6);
  mpz_add(Base->OldDenom, Temp7, Temp8);

  mpz_clear(Temp1);mpz_clear(Temp2);mpz_clear(Temp3);mpz_clear(Temp4);
  mpz_clear(Temp5);mpz_clear(Temp6);mpz_clear(Temp7);mpz_clear(Temp8);
}

#ifdef USE_BINARY_SPLITTING
void BinSplit(Point *Base, unsigned First, unsigned Last) {
   // Inclusive of First, Exclusive of Last
   Point Tip;

   switch (Last-First) {
      case 0:
         printf("BinSplit called with zero difference.\n");
         exit(0);
         break;
      case 1:
         BuildPoint(Base, First);
         break;
      default:
         InitPoint(&Tip);
         BinSplit(Base, First, (First+Last)/2);
         BinSplit(&Tip, (First+Last)/2, Last);
         MergePoints(Base, &Tip);
         DeInitPoint(&Tip);
         break;
   }
}
#endif

#ifdef USE_MATRIX_TOWER
Point Stack[100];
int StackNdx = -1;

void DumpStackSizes(void) {
   int x;
   printf("Stack sizes: ");
   for (x=0; x<=StackNdx; x++)
      printf("%d ", (int)mpz_sizeinbase(Stack[x].NewNum, 2));
   printf("\n");
}

Point Pop(void) {
   if (StackNdx < 0) {
      printf("Stack underflow.\n");
      exit(1);
   }
   return Stack[StackNdx--];
}

void Push(Point *Num) {
   if (StackNdx >=99) {
      printf("Stack overflow.\n");
      DumpStackSizes();
      exit(1);
   }
   Stack[++StackNdx] = *Num;
   /* We've pushed it, now re-init it so it wont point to the stack. */
   InitPoint(Num);
}

int ReadyToCollapse(void) {
   int Sz1, Sz2;

   if (StackNdx < 1)
      return 0;
   Sz1 = mpz_sizeinbase(Stack[StackNdx  ].NewNum, 2);
   Sz2 = mpz_sizeinbase(Stack[StackNdx-1].NewNum, 2);
   if (Sz1 >= Sz2)
      return 1;
   return 0;
}

void Collapse(void) {
   MergePoints(&Stack[StackNdx-1], &Stack[StackNdx]);
   DeInitPoint(&Stack[StackNdx]);
   StackNdx--;
}

int BuildTip(int NextTerm, int LastTerm) {
   // Inclusive NextTerm, Exclusive of LastTerm.
   // Do small chunks efficiently....
   int count, x;
   Point Tip, Term;

   InitPoint(&Tip);
   InitPoint(&Term);
   if (NextTerm >= LastTerm) {
      printf("Too far\n");
      exit(0);
   }
   // Aww, what the heck. Fake it.
   count = 8;
   if (count > (LastTerm-NextTerm))
      count = LastTerm - NextTerm;
   BuildPoint(&Tip, NextTerm);
   for (x=1; x<count; x++) {
      BuildPoint(&Term, NextTerm + x);
      MergePoints(&Tip, &Term);
   }
   Push(&Tip);
   DeInitPoint(&Tip);
   DeInitPoint(&Term);
   return count;
}
#endif

void ComputeCF(Point *Base, int NextTerm, int LastTerm,
               int a, int b, int c, int d) {
   /* starting matrix: (a*CF+b) / (c*CF+d) */
   /* Including NextTerm, excluding LastTerm */

   mpz_set_ui(Base->NewNum, a);
   mpz_set_ui(Base->OldNum, b);
   mpz_set_ui(Base->NewDenom, c);
   mpz_set_ui(Base->OldDenom, d);

#ifdef USE_BINARY_SPLITTING
   {
      Point Tip;
      InitPoint(&Tip);
      BinSplit(&Tip, NextTerm, LastTerm);
      MergePoints(Base, &Tip);
      DeInitPoint(&Tip);
   }
#endif
#ifdef USE_MATRIX_TOWER
   while (NextTerm < LastTerm) {
      if (ReadyToCollapse())
         Collapse();
      else
         NextTerm += BuildTip(NextTerm, LastTerm);
   }
   while (StackNdx > 0)
      Collapse();
   MergePoints(Base, &Stack[0]);
#endif
}

int main(int argc, char *argv[]) {
   int Digits = 0, Terms;
   Point pi;

   if (argc < 2) {
      printf("Need to know how many terms to process\n");
      exit(0);
   }
   Terms = atoi(argv[1]);
   InitPoint(&pi);
#if defined(NAT_LOG_PI)
   ComputeCF(&pi, 0, Terms+1, 0, 4, 1, 0); // 4/CF
   Digits = (int)(Terms/1.307) + 10;
#elif defined(NAT_LOG)
   ComputeCF(&pi, 0, Terms+1, 1, 0, 0, 1); // Identity
   Digits = (int)mpz_sizeinbase(pi.NewNum, 10) +
                 mpz_sizeinbase(pi.NewDenom, 10) + 10;
#elif defined(ARCCOTAN)
   ComputeCF(&pi, 0, Terms+1, 1, 0, 0, ARCCOTAN); // arctan(x)/x
   //ComputeCF(0, Terms+1, 1, 0, 0, 1); // Identity
   Digits = (int)(Terms*log10(4.0*ARCCOTAN_SQR)) + 10;
#elif defined(SQRT2)
   ComputeCF(&pi, 0, Terms+1, 1, 0, 0, 1); // Identity
   Digits = (int)mpz_sizeinbase(pi.NewNum, 10) +
                 mpz_sizeinbase(pi.NewDenom, 10) + 10;
#else
#error Unknown option.
#endif
   if (Digits == 0)
      Digits = Terms;
   printf("Terms %d Decimals %d\n", Terms, Digits);
#if 1
   // Print it out?
   {
      mpf_t N1, D1;
#define BITS_PER_DIGIT   3.32192809488736234787
      mpf_set_default_prec(Digits*BITS_PER_DIGIT);
      mpf_init(N1);
      mpf_init(D1);
      mpf_set_z(N1, pi.NewNum);
      mpf_set_z(D1, pi.NewDenom);
      mpf_div(N1, N1, D1);
      mpf_out_str(stdout, 10, Digits+10, N1); printf("\n");
      mpf_clear(N1);
      mpf_clear(D1);
   }
#endif
   DeInitPoint(&pi);
   return EXIT_SUCCESS;
}
