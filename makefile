# load module
#module load netcdf4
#module load gcc/4.8.0
# user library
OPTDIR		=	/home/cli/opt
# pooma library
POOMASUITE	=	gcc_openmp
POOMASRC	=	$(OPTDIR)/pooma241_$(POOMASUITE)/linux/src
POOMALIB	=	$(OPTDIR)/pooma241_$(POOMASUITE)/linux/lib
#jansson library
JANSSONLIB	=	$(OPTDIR)/jansson-2.6/lib
JANSSONINC	=	$(OPTDIR)/jansson-2.6/include
#netcdf library
NETCDFLIB	=	/opt/netcdf4/lib
NETCDFINC	=	/opt/netcdf4/include

CC		= g++
CFLAG	= -std=c++0x -O2 --openmp \
		  -ftemplate-depth-60 -fno-exceptions -Drestrict=__restrict__ -DNOPAssert -DNOCTAssert -funroll-loops
CLIB	= -L$(JANSSONLIB) -ljansson \
		  -L$(NETCDFLIB) -lnetcdf_c++ \
		  -L$(POOMALIB) -lpooma-gcc 
CINC	= -lm -I$(JANSSONINC) -I$(NETCDFINC) \
		  -I$(POOMALIB)/PoomaConfiguration-gcc -I$(POOMASRC) -I$(POOMALIB)\

MAIN 	= main
EXE		= run
SUPPORT = Boundary.h utils.h functions.h Stencil.h Variable.h setups.h configure.h\
		  Primitive.h \
		  NumericalScheme.h MicroPhysics.h

$(EXE): $(MAIN).cpp $(SUPPORT)
	$(CC) $(CFLAG) -o $(EXE) $(MAIN).cpp $(CINC) $(CLIB)

clean:
	@ rm -f $(EXE)
	@ rm -f *.o
.PHONY: clean
