
# build assorted tests for the Gnu MP code
#
# CC = cc -pg
CC = cc


EXES= polylog-bug unit-test 

MPLIB=../src/libanant.a
INC=../src

all: $(EXES)

polylog-bug.o: $(INC)/mp-binomial.h $(INC)/mp-complex.h \
               $(INC)/mp-misc.h $(INC)/mp-polylog.h 
unit-test.o: $(INC)/mp-zeta.h $(INC)/mp-binomial.h $(INC)/mp-complex.h \
             $(INC)/mp-consts.h $(INC)/mp-gamma.h $(INC)/mp-misc.h \
             $(INC)/mp-polylog.h $(INC)/mp-trig.h

polylog-bug:	polylog-bug.o $(MPLIB)

unit-test:	unit-test.o $(MPLIB)
	$(CC) -o unit-test $^ -lgmp -lgsl -lgslcblas -ldb -lm -lc

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

.c.o:
	$(CC) -c -g -O2 -Wall -I. -I../src $<

.o:
	$(CC) -g -o $* $^ $(MPLIB) -lgmp -ldb -lm
