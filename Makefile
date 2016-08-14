
IS_ICPC:= $(shell command -v icpc 2> /dev/null)

# STATIC_FLAG:= -static -static-intel

ifdef IS_ICPC
	CXX:=icpc
	CXXFLAGS:= -O3 -std=c++11 -xAVX -qopenmp -funroll-loops -ggdb
else
	CXX:=g++
	CXXFLAGS:= -O3 -std=c++11 -fopenmp -mtune=native -ggdb
endif

LD_INC_FLAGS:= -I./include 
LD_LIB_FLAGS:= -L./lib -lgfakluge  -lz -lhts #-lrodeo


SRC_DIR:=src
OBJ_DIR:=obj
LIB_DIR:=lib
INC_DIR:=include

EXE:=lasso

#$(EXE): $(SRC_DIR)/lasso.cpp $(LIB_DIR)/librodeo.a $(OBJ_DIR)/wrangler.o $(SRC_DIR)/vcfparse.hpp .pre-build
$(EXE): $(SRC_DIR)/lasso.cpp $(SRC_DIR)/vcfparse.hpp $(LIB_DIR)/libhts.a $(LIB_DIR)/libgfakluge.a .pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

$(LIB_DIR)/libhts.a: .pre-build
	+cd deps/htslib && $(MAKE) && cp libhts.a ../../lib/

$(LIB_DIR)/libvcflib.a: $(LIB_DIR)/libhts.a .pre-build
	+cd deps/vcflib && $(MAKE) && cp include/* ../../include/ && cp src/*.hpp ../../include/ && cp libvcflib.a ../../lib/

$(LIB_DIR)/libgfakluge.a: .pre-build
	+cd deps/gfakluge && $(MAKE) && cp gfakluge.hpp ../../include && cp libgfakluge.a ../../lib/

#$(LIB_DIR)/librodeo.a: $(OBJ_DIR)/rodeo.o $(OBJ_DIR)/wrangler.o .pre-build
#	ar -rs $@ $^

#$(OBJ_DIR)/rodeo.o: $(SRC_DIR)/rodeo.cpp $(SRC_DIR)/rodeo.hpp .pre-build
#	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

#$(OBJ_DIR)/wrangler.o: $(SRC_DIR)/wrangler.cpp $(SRC_DIR)/wrangler.hpp .pre-build $(LIB_DIR)/libgfakluge.a $(LIB_DIR)/libvcflib.a
#	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

.pre-build:
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	touch $@

clean:
	rm -rf $(OBJ_DIR)
	rm -rf $(LIB_DIR)
	+cd deps/vcflib && $(MAKE) clean
	+cd deps/gfakluge && $(MAKE) clean
	rm .pre-build

.PHONY: clean deps get-deps


