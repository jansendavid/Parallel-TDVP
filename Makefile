LIBRARY_DIR=$(ITENSOR)
include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk
INC+=-I${ITENSOR}
INC+=-I$(PWD)/include

LIBSPATH=-L$(ITENSOR)/lib
LIBSPATH+=$(LIBSLINK)



LIBS=-litensor 
LIBSG=-litensor-g 


#########################

CCFLAGS+=-I. $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS) -Wno-unused-variable -std=c++17 -O2 -std=gnu++1z
CCGFLAGS+=-I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS) 

LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS) -lboost_program_options -lboost_filesystem 
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS) -lboost_program_options -lboost_filesystem 
MPICOM=mpic++ -m64 -std=c++17 -fconcepts -fPIC


CPPFLAGS_EXTRA += -O2 -std=c++17



DB=-g
CXX=g++
ND=-DNDEBUG

##################################################



# tdmrg

execute_trottergates: src/execute_trottergates.cpp $(ITENSOR_LIBS) 
	$(CCCOM) $< -o bin/$@ $(CCFLAGS) $(INC) $(LIBFLAGS) 

# parallel tdvp
execute_p2tdvp: src/execute_p2tdvp.cpp $(ITENSOR_LIBS) 
	$(MPICOM) -std=c++17 $<  -o bin/$@  $(CCFLAGS) $(INC)  $(LIBSPATH) $(LIBFLAGS)  $(ND)



clean:
	rm bin/*

