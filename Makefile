CXX = g++
MARLEY_VERSION := $(shell cat ../marley-private/.VERSION)
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)
CXXFLAGS = -std=c++17 -Iinclude -I../marley-private/include -I/opt/homebrew/include \
  -DMARLEY_VERSION="\"$(MARLEY_VERSION)\"" -DUSE_ROOT $(ROOT_CFLAGS)
LDFLAGS = -L/opt/homebrew/lib -lgsl -lgslcblas $(ROOT_LIBS)
ALL_MARLEY_SRC = $(wildcard ../marley-private/src/*.cc)
MARLEY_SRC = $(filter-out ../marley-private/src/marley.cc ../marley-private/src/marsum.cc, $(ALL_MARLEY_SRC))

.PHONY: all clean run

ARGS :=

all: executables/omega executables/lambda

run: executables/lambda
	MARLEY=$(abspath ../marley-private) ./executables/lambda $(ARGS)

executables/omega: src/main_omega.cpp src/state_density.cpp src/preeq_common.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

executables/lambda: src/main_lambda.cpp src/transition_rates.cpp src/preeq_common.cpp src/state_density.cpp $(MARLEY_SRC)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f executables/omega executables/lambda