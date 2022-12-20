#HUGE_EDGE=1
#HUGE_VERTEX=1

ifdef HUGE_EDGE
INTE = -DHUGE_EDGE
endif

ifdef HUGE_VERTEX
INTV = -DHUGE_VERTEX
endif


CXX     = g++
CPPFLAGS = -O3 -std=c++11 -fopenmp -mavx -w $(INTE) $(INTV) -march=native -m64 -ftree-vectorize
DEBUG= -g -O0 -Wall  -std=c++11 -fopenmp -mavx -w 
LDFLAGS = -fopenmp -lpthread  
INC:=./include/
APP:=./app

target = pr bfs cc nibble sssp prdelta
debug = pr

all: $(target)

$(target): %: $(APP)/%.cpp $(DEP) 
	$(CXX) $(CPPFLAGS)  -o $(APP)/$@ $(DEP) -I $(INC) $<  $(LDFLAGS)

test: $(APP)/test.cpp
	$(CXX) $(CPPFLAGS)  -o $(APP)/$@  -I $(INC)  $<   $(LDFLAGS)

debug: $(APP)/$(debug).cpp $(DEP) 
	$(CXX) $(DEBUG)  -o $(APP)/$(debug) $(DEP) -I $(INC) $<  $(LDFLAGS)

clean:
	rm -f $(foreach bin, $(target), $(APP)/$(bin))