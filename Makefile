CXX = gcc

# Default configuration
CXXFLAGS_DEFAULT = -Ofast -Wall -Wno-float-equal -Wextra -std=c11 -march=native -fopenmp #-fopt-info-vec-missed 
LDFLAGS_DEFAULT = -lm -lquadmath

# Debug configuration
CXXFLAGS_DEBUG = -g -Wall -Wno-float-equal -Wextra
LDFLAGS_DEBUG = -lm -lquadmath

# Full debug configuration with sanitizers
CXXFLAGS_FULL_DEBUG = -g -Wall -Wno-float-equal -Wextra -std=c11 -fsanitize=address -S-fsanitize=undefined -fsanitize=leak
LDFLAGS_FULL_DEBUG = -lm -lquadmath

CXXFLAGS_PROFILE = -pg -Wall -Wno-float-equal -Wextra -std=c11
LDFLAGS_PROFILE = -lm  -lquadmath

TARGET = voronoi
SRC = helper.c force.c voronoi.c initial.c integrator.c logger.c thermo.c parser.c
HEADERS = jc_voronoi.h helper.h force.h voronoi.h initial.h integrator.h logger.h thermo.h parser.h

.PHONY: all default debug full_debug clean

all: default

default: $(TARGET)

debug: $(TARGET)_debug

full_debug: $(TARGET)_full_debug

profile: $(TARGET)_profile

$(TARGET): $(SRC) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS_DEFAULT) $(SRC) -o $(TARGET) $(LDFLAGS_DEFAULT)

$(TARGET)_debug: $(SRC) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS_DEBUG) $(SRC) -o $(TARGET) $(LDFLAGS_DEBUG)

$(TARGET)_full_debug: $(SRC) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS_FULL_DEBUG) $(SRC) -o $(TARGET) $(LDFLAGS_FULL_DEBUG)

$(TARGET)_profile: $(SRC) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS_PROFILE) $(SRC) -o $(TARGET) $(LDFLAGS_PROFILE)

clean:
	rm -f $(TARGET) $(TARGET)_debug $(TARGET)_full_debug