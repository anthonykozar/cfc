# c++ solaris/sunos/hpux compiler
#CCPP = CC 
#YACC = yacc
#LEX = lex

# c++ linux compiler (y.tab.h)
#CCPP = g++
#YACC = yacc
#LEX = lex

# cygwin compiler (cf_arith2.tab.h)
CCPP = g++
YACC = bison
LEX = flex

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
LIBS = libgmp.a libreadline.a libhistory.a -lfl -ltermcap # cygwin

cf_arith: $(OBJS) 
	$(CCPP) $(CFLAGS) -o $@ $(OBJS) $(LIBS) 

cf_arith.o: cf_arith.c cf_arith.h gmp.h cf_arith2.tab.h
	$(CCPP) $(CFLAGS) -c cf_arith.c

cf_arith2.tab.o: cf_arith2.tab.c cf_arith2.tab.h
	$(CCPP) $(CFLAGS) -c cf_arith2.tab.c

lex.yy.o: lex.yy.c cf_arith2.tab.h
	$(CCPP) $(CFLAGS) -c lex.yy.c

# dependencies

cf_arith2.tab.h cf_arith2.tab.c: cf_arith2.y gmp.h cf_arith.h
	$(YACC) $(YFLAGS) cf_arith2.y

lex.yy.c: cf_arith2.l gmp.h cf_arith.h
	$(LEX) cf_arith2.l

clean:
	/bin/rm -f *.o cf_arith2.tab.c cf_arith2.tab.h cf_arith2.output lex.yy.c cf_arith core cf_arith.exe cf_arith.exe.stackdump xout

