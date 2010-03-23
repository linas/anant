
# build Gnu MP code
#
# CC = cc -pg
CC = cc


EXES= db-merge db-prt unit-test 

MPLIB=libanant.a

all:  $(MPLIB) $(EXES)

MPOBJS= db-cache.o mp-binomial.o mp-cache.o mp-consts.o \
	mp-gamma.o mp-gkw.o mp-hyper.o mp-misc.o mp-polylog.o \
	mp-quest.o mp-trig.o mp-zeta.o

cache-fill:	cache-fill.o $(MPLIB)
db-merge:	db-merge.o $(MPLIB)
db-prt:	db-prt.o $(MPLIB)

# library objects
db-cache.o: db-cache.h
mp-binomial.o: mp-binomial.h mp-cache.h mp-complex.h mp-misc.h mp-trig.h
mp-cache.o: mp-cache.h mp-complex.h
mp-consts.o: mp-consts.h mp-binomial.h mp-complex.h mp-trig.h mp-zeta.h
mp-gamma.o: mp-gamma.h mp-binomial.h mp-complex.h mp-consts.h mp-misc.h mp-trig.h mp-zeta.h
mp-gkw.o: mp-gkw.h
mp-hyper.o: mp-hyper.h mp-complex.h mp-misc.h
mp-misc.o: mp-misc.h mp-complex.h
mp-polylog.o: mp-polylog.h mp-binomial.h mp-cache.h mp-complex.h mp-consts.h mp-gamma.h mp-misc.h mp-trig.h mp-zeta.h
mp-quest.o: mp-quest.h
mp-trig.o: mp-trig.h mp-binomial.h mp-cache.h mp-complex.h mp-misc.h
mp-zeta.o: mp-zeta.h db-cache.h mp-binomial.h mp-cache.h mp-complex.h mp-consts.h mp-trig.h

cache-fill.o: db-cache.h mp-zeta.h mp-misc.h
db-merge.o: db-cache.h mp-misc.h
unit-test.o: mp-zeta.h mp-binomial.h mp-complex.h mp-consts.h mp-gamma.h mp-misc.h mp-polylog.h mp-trig.h

unit-test:	unit-test.o $(MPLIB)
	$(CC) -o unit-test $^ -lgmp -lgsl -lgslcblas -ldb -lm -lc

$(MPLIB): $(MPOBJS)
	ar cru $(MPLIB) $^

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

.c.o:
	$(CC) -c -g -O2 -Wall $<

.o:
	$(CC) -g -o $* $^ $(MPLIB) -lgmp -ldb -lm