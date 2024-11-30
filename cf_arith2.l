%{
#include <setjmp.h>

#include "cf_arith.h"
#include "cf_arith2.tab.h"
%}

%%

-?[0-9]+ {                 /* integer with optional minus */
            mpz_init_set_str(yylval.mpzval, yytext, 10);
            return(NUMBER);
         }

scale    return(SCALE);    /* set precision e.g. scale = 50 */

output   return(OUTPUT);   /* set output format e.g. output = contfrac */

contfrac return(CONTFRAC); /* set output format e.g. output = contfrac */

decimal  return(DECIMAL);  /* set output format e.g. output = decimal */

rational return(RATIONAL); /* set output format e.g. output = rational */

timeout   return(TIMEOUT);   /* set timeout interrupt e.g. timeout = 100 */

cf       return(CF);       /* cf(xx) function */

sqrt     return(SQRT);     /* sqrt(xx) square root function */

atan     return(ATAN);     /* atan(xx) arc tangent function */

tan      return(TAN);      /* tan(xx)  tangent function */

log      return(LOG);      /* log(xx) logarithm function */

exp      return(EXP);      /* exp(xx) exponential function */


pi       return(PI);       /* pi continued fraction = [3 7 15 1 292...] */

e        return(E);        /* exp(1) continued fraction = [2 1 2 1 1 4...] */


x        return('x');      /* utilise dans {(1*x+3)/(1*x+4),pi} */

[A-Za-z][A-Za-z0-9]* {     /* variables */
         /* returns symbol's pointer */
         yylval.psymbolTable = findSymbol(yytext);
         return(VARIABLE);
                     }

[ \t]    ; /* ignore blanks */

"$"      return(0); /* end of input */

\n       | /* send to yacc */
.        { /* send to yacc */
            /* printf("**** %c\n", yytext[0]); */
            return(yytext[0]);
         }

%%

void scan_string(char *line) {
//#ifdef LINUX or cygwin
#if 1
   (void)yy_scan_string(line);
#endif
//#ifdef SOLARIS
#if 0
   char *xline = line;
   while (*xline) {
      xline++;
   }
   xline--;
   yyunput(EOF); // beginning at end
   while (xline >= line) {
      yyunput(*xline);
      xline--;
   }
#endif
}

void delete_string() {
//#ifdef LINUX
#if 0
   yy_delete_buffer(buffer_state);
#endif
//#ifdef cygwin
#if 1
   yy_delete_buffer(YY_CURRENT_BUFFER);
#endif
}

// (akozar) MOVED FROM cf_arith.c:
int yywrap(void) { return 1;}
