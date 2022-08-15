CXX     = g++
CPPFLAGS = -O3  -std=c++17 -fopenmp -mavx -march=native -m64 -ftree-vectorize -D_GLIBCXX_PARALLEL  -DINFO
DEBUG= -g -O0 -Wall  -std=c++17 -fopenmp -mavx -w -DDEBUG -D_GLIBCXX_PARALLEL -DINFO
LDFLAGS = -fopenmp  -lboost_timer -lboost_system -lboost_program_options # -ltbb -fopenmp-simd 
INC:=./include/
DEP:=#./include/partition.cpp #./include/moca.cpp ./include/propagate.tpp
APP:=./app

target = pr
#all: $(SOURCES) moca #pr add-spmv
debug = pr

all: $(target)

$(target): %: $(APP)/%.cpp
	$(CXX) $(CPPFLAGS) -o $(APP)/$@  -I $(INC)  $<  $(LDFLAGS)

test: $(APP)/test.cpp
	$(CXX) $(CPPFLAGS)  -o $(APP)/$@  -I $(INC)  $<   $(LDFLAGS)

debug: $(APP)/$(debug).cpp $(DEP) 
	$(CXX) $(DEBUG)  -o $(APP)/$(debug) $(DEP) -I $(INC) $<  $(LDFLAGS)

clean:
	rm -f $(APP)/$(target)

