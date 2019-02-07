include ${PETSC_DIR}/lib/petsc/conf/variables

CC=g++
VPATH=src
CFLAGS=-c -O3 -std=c++11 -I/home/adityak/Downloads/eigen-eigen-67e894c6cd8f -I./include -I${PETSC_DIR}/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include
OBJ := $(patsubst $(VPATH)/%.cpp,%.o,$(wildcard $(VPATH)/*.cpp))
CPP := $(patsubst $(VPATH)/%.cpp,%.cpp,$(wildcard $(VPATH)/*.cpp))

TARGET=FEM

$(TARGET): $(OBJ)
	$(CC) $^ -o $@ $(PETSC_LIB)
	rm *.o

$(OBJ): $(CPP)
	$(CC) $(CFLAGS) $^

clean:
	rm -f $(TARGET)
     
