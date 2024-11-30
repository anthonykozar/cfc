%{
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <setjmp.h>

#include "cf_arith.h"

// yylex() defined in lex.yy.c
extern int yylex();

extern int yyerror(char *); /* for cygwin/bison */

#define free(a) free((char *)a)

// global variables
struct chain *pcurrent;       // current pointer list
struct chain *pcurrentbegin;  // beginning of the current list
struct rationalcf *prational; // rational pointer
ContinuedFraction::enumOutput discriOutput; // decimal/continued fraction output

%}

%start xlist

%union {
   /* long int ival; */
   mpz_t mpzval;
   struct chain *pchain;
   struct rationalcf *prationalcf;
   ContinuedFraction *pvoid;
   struct symbolTable *psymbolTable;
}

/* terminal token type */
%token <psymbolTable> VARIABLE
%token <mpzval> NUMBER
%token CF SQRT ATAN TAN LOG EXP PI E SCALE OUTPUT CONTFRAC DECIMAL RATIONAL TIMEOUT

/* non terminal token type */
%type <pchain> list list_comma
%type <pvoid> cfexpr affectation
%type <prationalcf> fraction

%left '+' '-'
%left '*' '/'
%nonassoc UMINUS

%%

xlist      : setscale
           | setout
           | cfexpr
              {
                 if (discriOutput == ContinuedFraction::decimalOutput) {
                    // add one level for decimalOutput
                    // x = (1.x + 0)/(0.x + 1)
                    Matrix x(1, 0, 0, 1);
                    ContinuedFraction *cf1 = new ContinuedFraction(x,
                                                    *$1, ContinuedFraction::decimalOutput);
                    (*cf1).x0Andx1Init();
                    (*cf1).evalMain();
                 } else if (discriOutput == ContinuedFraction::cfOutput) {
                    (*$1).evalMain();
                 } else if (discriOutput == ContinuedFraction::rationalOutput) {
                    // add one level for rationalOutput
                    // x = (1.x + 0)/(0.x + 1)
                    Matrix x(1, 0, 0, 1);
                    ContinuedFraction *cf1 = new ContinuedFraction(x, *$1);
                    //(*cf1).x0Andx1Init();
                    (*cf1).evalMain();
                 }
              }
           | affectation
              {
              }
           | error
              {
                 //yyerrok;
              }
           ;

setscale   : SCALE '=' NUMBER
              {
                 if (mpz_cmp_ui($3, 0) >= 0) {
                    mpz_set(scale, $3);
                    printf("ok\n");
                 } else {
                    printf("ko\n");
                 }
              }
           ;

setout     : OUTPUT '=' CONTFRAC
              {
                 discriOutput = ContinuedFraction::cfOutput;
                 printf("ok\n");
              }
           | OUTPUT '=' DECIMAL
              {
                 discriOutput = ContinuedFraction::decimalOutput;
                 printf("ok\n");
              }
           | OUTPUT '=' RATIONAL
              {
                 discriOutput = ContinuedFraction::rationalOutput;
                 printf("ok\n");
              }
           ;

cfexpr     : '(' cfexpr ')'
              {
                 $$ = $2;
              }
           | cfexpr '+' cfexpr
              {
                 ContinuedFraction *cf1 = new ContinuedFraction();
                 *cf1 = (*cf1).plus(*$1, *$3, ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | cfexpr '-' cfexpr
              {
                 ContinuedFraction *cf1 = new ContinuedFraction();
                 *cf1 = (*cf1).substract(*$1, *$3, ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | cfexpr '*' cfexpr
              {
                 ContinuedFraction *cf1 = new ContinuedFraction();
                 *cf1 = (*cf1).multiply(*$1, *$3, ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | cfexpr '/' cfexpr
              {
                 ContinuedFraction *cf1 = new ContinuedFraction();
                 *cf1 = (*cf1).divide(*$1, *$3, ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | '-' cfexpr %prec UMINUS
              {
                 // x = (-1.x + 0)/(0.x + 1)
                 Matrix x(-1, 0, 0, 1);
                 ContinuedFraction *cf1 = new ContinuedFraction(x,
                                                 *$2, ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | '[' list_comma ']'
              {
                 ContinuedFraction *cf1 = new ContinuedFraction($2,
                                                 ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | CF '(' fraction ')'   /* for instance : cf(xx/yy) */
              {
                 ContinuedFraction *cf1 = new ContinuedFraction(
                                                 $3->num, $3->den,
                                                 ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | SQRT '(' fraction ')' /* for instance : sqrt(xx/yy) */
              {
                 if (((mpz_cmp_ui($3->num, 0) < 0) && (mpz_cmp_ui($3->den, 0) > 0)) ||
                     ((mpz_cmp_ui($3->num, 0) > 0) && (mpz_cmp_ui($3->den, 0) < 0))) {
                    printf("error : negative argument for sqrt()\n");
                    longjmp(env, 1);
                 } else {
                    ContinuedFraction *cf1 = new ContinuedFraction();
                    *cf1 = (*cf1).squareRootRat($3->num, $3->den,
                                                ContinuedFraction::cfOutput);
                    (*cf1).x0Andx1Init();
                    $$ = cf1;
                 }
              }
           | SQRT '(' cfexpr ')'   /* for instance : sqrt([xxx]) */
              {
                 // A MODIFIER !!! tester le cas ou l'argument est negatif !!!
                 ContinuedFraction *cf1 = new ContinuedFraction();
                 *cf1 = (*cf1).squareRootOp(*$3, ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | EXP '(' fraction ')' /* rational exponentiation : exp(xx/yy) */
              {
                 ContinuedFraction *cf1 = new ContinuedFraction(expRat,
                                                 $3->num, $3->den,
                                                 ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | LOG '(' fraction ')' /* rational logarithm : log(xx/yy) */
              {
                 if (((mpz_cmp_ui($3->num, 0) < 0) && (mpz_cmp_ui($3->den, 0) > 0)) ||
                     ((mpz_cmp_ui($3->num, 0) > 0) && (mpz_cmp_ui($3->den, 0) < 0))) {
                    printf("error : negative argument for log()\n");
                    longjmp(env, 1);
                 } else if (mpz_cmp($3->num, $3->den) < 0) {
                    // if n/d < 1 then log(n/d) = - log (d/n) 
                    // add one level for -log(1/x)
                    ContinuedFraction *cf1 = new ContinuedFraction(logRat,
                                                    $3->den, $3->num,
                                                    ContinuedFraction::cfOutput);
                    (*cf1).x0Andx1Init();
                    // -x = (-1.x + 0)/(0.x + 1)
                    Matrix x(-1, 0, 0, 1);
                    ContinuedFraction *cf2 = new ContinuedFraction(x,
                                                    *cf1, ContinuedFraction::cfOutput);
                    (*cf2).x0Andx1Init();
                    $$ = cf2;
                 } else {
                    ContinuedFraction *cf1 = new ContinuedFraction(logRat,
                                                    $3->num, $3->den,
                                                    ContinuedFraction::cfOutput);
                    (*cf1).x0Andx1Init();
                    $$ = cf1;
                 }
              }
           | ATAN '(' fraction ')' /* rational arc tangent : atan(xx/yy) */
              {
                 ContinuedFraction *cf1 = new ContinuedFraction(atanRat,
                                                 $3->num, $3->den,
                                                 ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | TAN '(' fraction ')' /* rational tangent : tan(xx/yy) */
              {
                 ContinuedFraction *cf1 = new ContinuedFraction(tanRat,
                                                 $3->num, $3->den,
                                                 ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | PI
              {
                 ContinuedFraction *cf1 = new ContinuedFraction(&pi,
                                                 ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | E
              {
                 ContinuedFraction *cf1 = new ContinuedFraction(&e,
                                                 ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | '{' '(' NUMBER '*' 'x' '+' NUMBER ')' '/' '(' NUMBER '*' 'x' '+' NUMBER ')' ',' cfexpr '}'
              {
                 Matrix x($3, $7, $11, $15);
                 ContinuedFraction *cf1 = new ContinuedFraction(x, *$18,
                                                 ContinuedFraction::cfOutput);
                 (*cf1).x0Andx1Init();
                 $$ = cf1;
              }
           | VARIABLE
              {
                 ContinuedFraction *cf1 = $1->value;
                 $$ = cf1;
              }
           ;

fraction   : NUMBER '/' NUMBER     /* for instance : 355/113 */
              {
                 // A MODIFIER !!! tester le cas num<0 et den<0
                 prational = (struct rationalcf *)malloc(sizeof(struct rationalcf));
                 mpz_init(prational->num);
                 mpz_init(prational->den);
                 mpz_set(prational->num, $1);
                 mpz_set(prational->den, $3);
                 if (mpz_cmp_ui(prational->den, 0) == 0) {
                    printf("error : null denominator\n");
                    longjmp(env, 1);
                 } else {
                    $$ = prational;
                 }
              }
           | NUMBER                /* for instance : 163 */
              {
                 prational = (struct rationalcf *)malloc(sizeof(struct rationalcf));
                 mpz_init(prational->num);
                 mpz_init(prational->den);
                 mpz_set(prational->num, $1);
                 mpz_set_ui(prational->den, 1);
                 $$ = prational;
              }
           | NUMBER '.' NUMBER     /* for instance : 3.141592 */
              {
                 prational = (struct rationalcf *)malloc(sizeof(struct rationalcf));
                 mpz_init(prational->num);
                 mpz_init(prational->den);
                 // not finished
                 $$ = prational;
              }
           ;

list_comma : list ',' list
              {
                 struct chain *pchain = $1;
                 while (pchain->next != NULL) {
                    pchain = pchain->next;
                 }
                 pchain->next = $3; /* concatenation */
                 pchain = $3;
                 while (pchain->next != NULL) {
                    pchain = pchain->next;
                 }
                 pchain->next = $3; /* loop */
                 $$ = $1;
              }
           | ',' list
              {
                 struct chain *pchain = $2;
                 while (pchain->next != NULL) {
                    pchain = pchain->next;
                 }
                 pchain->next = $2; /* loop */
                 $$ = $2;
              }
           | list
              {
                 struct chain *pchain = $1;
                 $$ = $1;
              }
           ;

list       : NUMBER
              {
#ifdef DEBUG
                 printf("*** list : NUMBER\n");
                 mpz_out_str(stdout, 10, $1); printf("\n");
#endif
                 pcurrent = (struct chain *)malloc(sizeof(struct chain));
                 mpz_init(pcurrent->value);
                 mpz_set(pcurrent->value, $1);
                 pcurrent->next = NULL;
                 pcurrentbegin = pcurrent;
                 $$ = pcurrentbegin;
              }
           | list NUMBER
              {
#ifdef DEBUG
                 printf("*** list : list | NUMBER\n");
                 mpz_out_str(stdout, 10, $2); printf("\n");
#endif
                 pcurrent->next = (struct chain *)malloc(sizeof(struct chain));
                 pcurrent = pcurrent->next;
                 mpz_init(pcurrent->value);
                 mpz_set(pcurrent->value, $2);
                 pcurrent->next = NULL;
                 $$ = pcurrentbegin;
              }
           ;

affectation : VARIABLE '=' cfexpr
              {
                 printf("*** affectation\n");
                 // TODO : not a TRUE affectation : make a new and then affect...
                 $1->value = $3;
              }
           ;

%% /* start of programs */

int yyerror(char *s) {
   fprintf(stderr, "*** %s *** yytname = <%s>\n", s, yytname[YYTRANSLATE(yychar)]);
   return 0;
}

#ifdef DEBUG
void xprint(struct chain *pchain) {
   int i=0;

   printf("*** xprint()\n");
   printf("*** [ ");
   while (pchain != NULL) {
      mpz_out_str(stdout, 10, pchain->value); printf(" ");
      pchain = pchain->next;
      i++; if (i==50) { printf("...] ***\n"); return; }
   }
   printf("] ***\n");
}
#endif
