PYTHONINCLUDE = $(shell python -c "import sys; print sys.prefix+'/include/python'+sys.version[:3]")
PYTHONLIB = $(shell python -c "import sys; print sys.prefix+'/lib/python'+sys.version[:3]")
PYTHONSITEPKGS = $(PYTHONLIB)/site-packages
DESTDIR = $(PYTHONSITEPKGS)/orange
OS = $(shell uname)

COMPILER = g++
CCOMPILER = gcc
LINKER = g++

LIBRARY_FILE = _orngCRS.so
DEPS = clustind_wrap.o twins.o pam.o fanny.o logreg.o lsq.o interface.o svm.o nb.o kikuchi.o pr_statistics.o models.o

ifeq ($(OS), Darwin)
  COMPILEOPTIONS = -arch i386 -arch ppc -fPIC -fpermissive -fno-common -w -O3 -I$(PYTHONINCLUDE)
  LINKOPTIONS = -arch i386 -arch ppc -dynamiclib -headerpad_max_install_names -undefined dynamic_lookup -lstdc++
else
  COMPILEOPTIONS = -fPIC -fpermissive -w -O3 -I$(PYTHONINCLUDE)
  LINKOPTIONS = -shared -lstdc++
endif

all: $(LIBRARY_FILE)

$(LIBRARY_FILE): $(DEPS)
	$(LINKER) $(LINKOPTIONS) $(DEPS) -o $(LIBRARY_FILE)
ifeq ($(OS), Darwin)
	install_name_tool -id $(DESTDIR)/$(LIBRARY_FILE) $(LIBRARY_FILE)
endif

clustind_wrap.o: Clustind_wrap.cxx
	$(COMPILER) -c $(COMPILEOPTIONS) Clustind_wrap.cxx -o clustind_wrap.o

Clustind_wrap.cxx: clustind.i
	swig -python clustind.i

%.o: %.cpp
	$(COMPILER) $(COMPILEOPTIONS) -c $< -o $@

%.o: %.c
	$(CCOMPILER) $(COMPILEOPTIONS) -c $< -o $@

clean:
	rm -f *.o $(LIBRARY_FILE)

.PHONY: all clean
