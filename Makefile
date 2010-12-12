

VERSION=0.2.0
VERDIR=anant-${VERSION}

dist:
	rm -rf ${VERDIR}
	mkdir ${VERDIR}
	mkdir ${VERDIR}/src
	cp LICENSE ${VERDIR}
	cp Makefile ${VERDIR}
	cp README ${VERDIR}
	cp src/*.c ${VERDIR}/src
	cp src/*.h ${VERDIR}/src 
	cp src/Makefile ${VERDIR}/src
	tar -zcvf ${VERDIR}.tar.gz ${VERDIR}
