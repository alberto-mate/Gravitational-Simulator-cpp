# Compile all the files in the directory using the make command


CC = g++
CFLAGS = -Wall -Wextra -fopenmp
LDFLAGS = `pkg-config --cflags --libs opencv4`

SRC_DIR = .
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:.cpp=.o)

all: $(OBJS)

%.o: %.cpp
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(OBJS)