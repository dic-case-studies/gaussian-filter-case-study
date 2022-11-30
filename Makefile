CXXFLAGS=-std=c++14 -Wall -Wextra -pedantic -I include -O3 `pkg-config --libs cfitsio`

processor := $(shell uname -m)
ifeq ($(processor),$(filter $(processor),aarch64 arm64))
    ARCH_C_FLAGS += -march=armv8-a+fp+simd+crc 
else ifeq ($(processor),$(filter $(processor),i386 x86_64))
    ARCH_C_FLAGS += -march=native 
endif

DEBUGFLAGS=-fsanitize=address -g

CXX=g++

all: dir build/main 

build/%: src/%.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(ARCH_C_FLAGS) $(OPT)

scratch: scratch.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(ARCH_C_FLAGS) $(OPT)

clean:
	rm -rf build/* *app scratch

dir:
	mkdir -p build

.PHONY: all clean

