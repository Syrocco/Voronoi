CXX = gcc

# Default configuration
CXXFLAGS_DEFAULT = -Ofast -Wall -Wno-float-equal -Wextra -std=c11 -march=native
LDFLAGS_DEFAULT = -lm

# Debug configuration
CXXFLAGS_DEBUG = -g -Wall -Wno-float-equal -Wextra -std=c11
LDFLAGS_DEBUG = -lm

# Full debug configuration with sanitizers
CXXFLAGS_FULL_DEBUG = -g -Wall -Wno-float-equal -Wextra -std=c11 -fsanitize=address -fsanitize=undefined -fsanitize=leak
LDFLAGS_FULL_DEBUG = -lm 

TARGET = voronoi
SRC = voronoi.c helper.c force.c
HEADERS = helper.h force.h voronoi.h

.PHONY: all default debug full_debug clean

all: default

default: $(TARGET)

debug: $(TARGET)_debug

full_debug: $(TARGET)_full_debug

$(TARGET): $(SRC) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS_DEFAULT) $(SRC) -o $(TARGET) $(LDFLAGS_DEFAULT)

$(TARGET)_debug: $(SRC) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS_DEBUG) $(SRC) -o $(TARGET) $(LDFLAGS_DEBUG)

$(TARGET)_full_debug: $(SRC) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS_FULL_DEBUG) $(SRC) -o $(TARGET) $(LDFLAGS_FULL_DEBUG)

clean:
	rm -f $(TARGET) $(TARGET)_debug $(TARGET)_full_debug