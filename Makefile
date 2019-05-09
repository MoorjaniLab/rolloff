HOMEL=$(PWD)
TOP=bin
## binary for install  
DEBUG_OPTIONS= -g -pg ## -fprofile-arcs -ftest-coverage 
BIN=$(HOMEL)/bin
PERLDIR=$(HOMEL)/perlsrc

# "make clean" to clean up extra files in this directory
# "make clobber" to clobber all files and subdirectories except source files
#    so as to enable recompiling from scratch.

NLIB=$(HOMEL)/smartlib/nicklib.a
IDIR=$(HOMEL)/smartinclude
IDIR2=$(HOMEL)/mcmcsrc
VPATH=.:nicksrc:mcmcsrc:perlsrc
BLAS = blas
# may need to change to BLAS = blas (depends on blas/lapack installation)

###CC=/util/bin/gcc 
FF=gfortran     
##FF=g77

CFLAGS= -c -g -I$(IDIR) -I$(IDIR2)  -Wimplicit ## -fprofile-arcs -ftest-coverage 
OBJ=strsubs.o sortit.o vsubs.o statsubs.o linsubs.o getpars.o xsearch.o gauss.o	gds.o
TWTAB=\"$(HOME)/tables/twtable\"


statsubs.o:     nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -DTWTAB=$(TWTAB) -o statsubs.o nicksrc/statsubs.c

M1=rolloff
M1O=rolloff.o  admutils.o  ldsubs.o  mcio.o  regsubs.o  egsubs.o

$(M1): nicklib $(M1O)
	gcc -O -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M1) $(M1O) $(NLIB) -Wimplicit


libnick.a:	dirs tables  $(OBJ)
	ar -r libnick.a $(OBJ)

nicklib:	libnick.a 
	cp libnick.a  $(NLIB)

tables:    
	echo "tables made"  > tables
	cp twtable  $(HOMEL)/smarttables

dirs:	
	mkdir -p  $(HOMEL)/smartlib
	mkdir -p  $(HOMEL)/smarttables
	mkdir -p  $(HOMEL)/smartinclude
	mkdir -p  $(BIN)
	cp  *.h  $(IDIR)
	cp  nicksrc/*.h  $(IDIR)

clean: 
	rm -f *.o 
	rm -f *junk*
	rm -f core
	rm -f libnick.a
	rm -f $(PROGS)
	rm -f nicksrc/*.o

clobber: clean rmdirs 


rmdirs: 
	rm -rf $(HOMEL)/smartlib 
	rm -rf $(HOMEL)/smarttables 
	rm -rf $(HOMEL)/smartinclude 

