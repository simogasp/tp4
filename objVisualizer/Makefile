CC = g++
CXXFLAGS = -Wall -O3 -DHAVE_STD_UNORDERED_MAP -DHAVE_STD_FUNCTIONAL -std=c++0x
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LIBS := -framework OpenGL -lglut -lm
else
LIBS = -lGL -lglut -lm -lGLU
endif


all: main 

# solutions
main: main.o ObjModel.o core.o
	$(CC)  $(CXXFLAGS) -o $@ $^ $(LIBS)

solution: main.o ObjModelsolution.o core.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o main solution
