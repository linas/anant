

VERSION=0.2.1
VERDIR=anant-${VERSION}

dist:
	rm -rf ${VERDIR}
	mkdir ${VERDIR}
	mkdir ${VERDIR}/src
	mkdir ${VERDIR}/tests
	cp LICENSE ${VERDIR}
	cp Makefile ${VERDIR}
	cp README ${VERDIR}
	cp src/*.c ${VERDIR}/src
	cp src/*.h ${VERDIR}/src 
	cp src/Makefile ${VERDIR}/src
	cp tests/*.c ${VERDIR}/tests
	cp tests/*.h ${VERDIR}/tests 
	cp tests/Makefile ${VERDIR}/tests
	tar -zcvf ${VERDIR}.tar.gz ${VERDIR}
