CXX = gcc
CXXFLAGS = -g -Wall -Wno-float-equal -Wextra -std=c11 
LDFLAGS = -lm

TARGET = voronoi
SRC = voronoi.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET)