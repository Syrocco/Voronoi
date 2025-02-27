CXX = gcc
CXXFLAGS = -O3 -Wall -Wextra -std=c11 
LDFLAGS = -lm

TARGET = voronoi
SRC = voronoi.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET)