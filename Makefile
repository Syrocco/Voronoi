CXX = gcc
CXXFLAGS = -g -Wall -Wno-float-equal -Wextra -std=c11 #-fsanitize=address -fsanitize=undefined -fsanitize=leak 
LDFLAGS = -lm

TARGET = voronoi
SRC = voronoi.c helper.c force.c

all: $(TARGET)

$(TARGET): $(SRC) 
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET)
