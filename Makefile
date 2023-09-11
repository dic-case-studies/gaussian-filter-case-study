CXXFLAGS=-std=c++14 -Wall -Wextra -pedantic -I include `pkg-config --libs --cflags cfitsio` -DASSERT

processor := $(shell uname -m)
ifeq ($(processor),$(filter $(processor),aarch64 arm64))
    ARCH_C_FLAGS += -march=armv8-a+fp+simd+crc -D__SSE__
else ifeq ($(processor),$(filter $(processor),i386 x86_64))
    ARCH_C_FLAGS += -march=native 
endif

OPTIMISATIONS=-O3
# OPTIMISATIONS=-fsanitize=address -g

CXX=g++

all: dir build/main 

src/%: src/%.cpp dir
	$(CXX) $(CXXFLAGS) $(ARCH_C_FLAGS) $(OPTIMISATIONS) $(OPT) -o build/$@ $<

clean:
	rm -rf build/*

dir:
	mkdir -p build/src

.PHONY: all clean

