# densest lattice packings
all:
	make densest_lattice_packings

CC = gcc
CXX = g++
CXXFLAGS =  -Wp,-MMD,$*.dep \
            -Wall -Wextra -pedantic \
            -fPIC \
            -Wl,-no_pie \
            -g

# adjust these paths to point to your soplex folder!
#SOPLEX_INC   = -I../soplex-1.7.0/src
#SOPLEX_LIB = -L../soplex-1.7.0/lib -lsoplex -lgmp -lz 

#PM_INC := $(shell /store/borel/lorenz/cluster/prefix/bin/polymake-config --includes)
#PM_CFLAGS := -fPIC -DPOLYMAKE_DEBUG=0
#PM_LIBS := $(shell /store/borel/lorenz/cluster/prefix/bin/polymake-config --libs)
#PM_LDFLAGS := $(shell /store/borel/lorenz/cluster/prefix/bin/polymake-config --ldflags)

CFLAGS     = -O3 -Wall -fno-strict-aliasing -pg

OBJ = main.o error.o helper.o aabox.o polytope.o densest_lattice_packing.o numerics.o linear_programming.o admissibleLattice.o LPSolver.o Matrix.o Vector.o Polynomial3.o Minima.o

OBJ_RAND = randomize.o error.o helper.o aabox.o polytope.o densest_lattice_packing.o numerics.o linear_programming.o admissibleLattice.o LPSolver.o Matrix.o Vector.o Polynomial3.o Minima.o

OBJ_ROOT_CASE = main_wip_test.o linear_factor.o Polynomial3.o helper.o

LDFLAGS = -L../pari-2.7.4/Odarwin-i386
LIB = -ldl -lm -lc -L/usr/local/lib -lpari 
#$(SOPLEX_LIB)
INC =  -I. -I../pari-2.7.4/src/headers ${PM_INC}
#$(SOPLEX_INC)

densest_lattice_packings: $(OBJ)
	$(CXX) -o $@ $(OBJ) ${LDFLAGS} ${PM_LDFLAGS} ${LIB} ${PM_LIBS} -Wl,--warn-unresolved-symbols ${LDFLAGS} -pg

root_case_test: $(OBJ_ROOT_CASE)
	$(CXX) -o $@ $(OBJ_ROOT_CASE) ${LDFLAGS} -Wl ${LDFLAGS} -pg

main.o: main.cpp main.h
	$(CXX) -o main.o -c -I. $(INC) ${PM_INC} ${PM_CFLAGS} ${CFLAGS} main.cpp

main_wip_test.o: main_wip_test.cpp main_wip_test.h helper.o
	$(CXX) -o main_wip_test.o -c -I. $(INC)  ${CFLAGS} main_wip_test.cpp

linear_factor.o: linear_factor.cpp linear_factor.h
	$(CXX) -o linear_factor.o -c -I. $(INC)  ${CFLAGS} linear_factor.cpp

numerics.o: numerics.c numerics.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC numerics.c

admissibleLattice.o: admissibleLattice.cpp admissibleLattice.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC admissibleLattice.cpp

error.o: error.cpp error.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC error.cpp

helper.o: helper.cpp helper.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC helper.cpp

aabox.o: aabox.cpp aabox.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC aabox.cpp

polytope.o: polytope.cpp polytope.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC polytope.cpp

densest_lattice_packing.o: densest_lattice_packing.cpp densest_lattice_packing.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC densest_lattice_packing.cpp

linear_programming.o: linear_programming.cpp linear_programming.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC linear_programming.cpp
    
LPSolver.o: LPSolver.cpp LPSolver.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC LPSolver.cpp
    
Matrix.o: Matrix.cpp Matrix.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC Matrix.cpp
    
Vector.o: Vector.cpp Vector.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC Vector.cpp
    
Polynomial3.o: Polynomial3.cpp Polynomial3.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC Polynomial3.cpp
    
Minima.o: Minima.cpp Minima.h
	$(CXX) -c $(INC) $(CFLAGS) -fPIC Minima.cpp

dlplpfileconverter: dlplpfileconverter.o
	$(CXX) -o $@ dlplpfileconverter.o ${PM_LDFLAGS} ${PM_LIBS} -lc 

dlplpfileconverter.o: dlplpfileconverter.cpp main.h polytope.h polytope.cpp
	$(CXX) -c main.h polytope.h $(CFLAGS) -fPIC dlplpfileconverter.cpp
    
randomized_packings: $(OBJ) randomize.o
	$(CXX) -o $@ $(OBJ_RAND) ${PM_LDFLAGS} ${LIB} ${PM_LIBS} -pg

randomize.o: randomize.cpp
	$(CXX) -c $(INC) $(CFLAGS) -fPIC randomize.cpp


clean:
	-rm *.o

