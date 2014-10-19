CC = g++
CCOPTS = -Wall -O3

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LIBS := -framework OpenGL -lglut -lm
else
LIBS = -lGL -lglut -lm -lGLU
endif


all: glutobj 

# solutions
glutobj: glutobj.o ObjModel.o
	$(CC) $(CCOPTS) -o $@ $^ $(LIBS)

glutobjsolution: glutobjsolution.o ObjModelsolution.o
	$(CC) $(CCOPTS) -o $@ $^ $(LIBS)
	
clean:
	rm -f *.o glutobj glutobjsolution

