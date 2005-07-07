#INSTALL NOTES: Before running this makefile, change the INCLUDE_DIR and 
#PYTHONLIB variables to for your local environment

INCLUDE_DIR = /usr/include/python2.4
PYTHONLIB = /usr/lib/python2.4/config/libpython2.4.dll.a

CFLAGS = -fPIC


default:
	make all

clean:
	rm -f *.o *.so

all: _orngCRS

clustind_wrap.cxx: 
	swig -python clustind.i

clustind_wrap.o: 
	g++ -c $(CFLAGS) clustind_wrap.cxx -o clustind_wrap.o -I$(INCLUDE_DIR)

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

models.o: models.cpp
	g++ -c $(CFLAGS) models.cpp -o models.o

kikuchi.o: kikuchi.cpp
	g++ -c $(CFLAGS) kikuchi.cpp -o kikuchi.o

nb.o: nb.cpp
	g++ -c $(CFLAGS) nb.cpp -o nb.o

pr_statistics.o: pr_statistics.cpp
	g++ -c $(CFLAGS) pr_statistics.cpp -o pr_statistics.o


interface.o: interface.cpp
	g++ -c $(CFLAGS) interface.cpp -o interface.o -I$(INCLUDE_DIR)

_orngCRS: clustind_wrap.o twins.o pam.o fanny.o logreg.o lsq.o interface.o svm.o nb.o kikuchi.o pr_statistics.o models.o
	g++ -shared clustind_wrap.o twins.o pam.o fanny.o logreg.o lsq.o interface.o svm.o nb.o kikuchi.o pr_statistics.o models.o $(PYTHONLIB) -lm -o _orngCRS.so
	cp -f _orngCRS.so ../bin

