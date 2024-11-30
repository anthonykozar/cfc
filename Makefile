# Mac OS X changes 2018-06-17 by Anthony Kozar:
#   - added INCL option for the location of gmp.h
#   - removed gmp.h as a dependency in the build rules
#   - renamed local gmp.h to old-gmp.h
#   - moved definition of yywrap() from cf_arith.c to cf_arith.l
#   - TO DO: change extensions for cf_arith.c, cf_arith2.tab.c, lex.yy.c to .cpp ?

# c++ solaris/sunos/hpux compiler
#CCPP = CC 
#YACC = yacc
#LEX = lex

# c++ linux compiler (y.tab.h)
#CCPP = g++
#YACC = yacc
#LEX = lex

# cygwin compiler (cf_arith2.tab.h)
#CCPP = g++
#YACC = bison
#LEX = flex

# Mac OS X compiler (cf_arith2.tab.h)
CCPP = g++
YACC = bison
LEX = flex
INCL = -I/usr/local/includes

# -g : info pour debug
###/usr/nlang/dbx cf_arith core
###(dbx) where
###(dbx) quit
#CFLAGS = -DDEBUG # -DXDEBUG
CFLAGS = -g

# -d : force creation of cf_arith2.tab.h
# -v : verbose option to obtain cf_arith2.output 
# -t : debug option
YFLAGS = -d -v -t

OBJS = cf_arith.o cf_arith2.tab.o lex.yy.o

# liby.a contains main() and yyerror()
#LIBS = libgmp.a libreadline.a libhistory.a -ll -ltermcap # -ly # sunos
#LIBS = libgmp.a -lfl       # linux
#LIBS = libgmp.a libreadline.a libhistory.a -lfl -ltermcap # cygwin
LIBS = -lgmp -lreadline -ltermcap # Mac OS X

cf_arith: $(OBJS) 
	$(CCPP) $(CFLAGS) -o $@ $(OBJS) $(LIBS) 

cf_arith.o: cf_arith.c cf_arith.h cf_arith2.tab.h # gmp.h
	$(CCPP) $(CFLAGS) $(INCL) -c cf_arith.c

cf_arith2.tab.o: cf_arith2.tab.c cf_arith2.tab.h
	$(CCPP) $(CFLAGS) $(INCL) -c cf_arith2.tab.c

lex.yy.o: lex.yy.c cf_arith2.tab.h
	$(CCPP) $(CFLAGS) $(INCL) -c lex.yy.c

# dependencies

cf_arith2.tab.h cf_arith2.tab.c: cf_arith2.y cf_arith.h # gmp.h
	$(YACC) $(YFLAGS) cf_arith2.y

lex.yy.c: cf_arith2.l cf_arith.h # gmp.h
	$(LEX) cf_arith2.l

clean:
	/bin/rm -f *.o cf_arith2.tab.c cf_arith2.tab.h cf_arith2.output lex.yy.c cf_arith core cf_arith.exe cf_arith.exe.stackdump xout

