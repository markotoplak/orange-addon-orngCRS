#INSTALL NOTES: Before running this makefile, change the INCLUDE_DIR and 
#PYTHONLIB variables to for your local environment

INCLUDE_DIR = $(HOME)/include/python2.2
PYTHONLIB = $(HOME)/lib/libpython2.2.so

CFLAGS = -fPIC


default:
	make all

clean:
	rm -f *.o *.so

all: _orngCRS

clustind_wrap.c: 
	swig -python clustind.i

clustind_wrap.o: 
	gcc -c $(CFLAGS) clustind_wrap.c -o clustind_wrap.o -I$(INCLUDE_DIR)

twins.o: twins.c
	gcc -c $(CFLAGS) twins.c -o twins.o

pam.o: pam.c
	gcc -c $(CFLAGS) pam.c -o pam.o 

fanny.o: fanny.c
	gcc -c $(CFLAGS) fanny.c -o fanny.o 

logreg.o: logreg.cpp
	g++ -c $(CFLAGS) logreg.cpp -o logreg.o

lsq.o: lsq.cpp
	g++ -c $(CFLAGS) lsq.cpp -o lsq.o

svm.o: svm.cpp
	g++ -c $(CFLAGS) svm.cpp -o svm.o

interface.o: interface.cpp
	g++ -c $(CFLAGS) interface.cpp -o interface.o -I$(INCLUDE_DIR)

_orngCRS: clustind_wrap.o twins.o pam.o fanny.o logreg.o lsq.o interface.o svm.o
	g++ -shared clustind_wrap.o twins.o pam.o fanny.o logreg.o lsq.o interface.o svm.o -lm -L$(PYTHONLIB) -o _orngCRS.so
	cp -f _orngCRS.so ../bin

