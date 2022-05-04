
#Compiler and project folder flags
CC = g++
FLAGS = -std=c++17

SRC = src
OBJ = obj
RES = results
BIN = bin/output

#List of project subdiretories
DIRS = $(filter-out $(SRC)/main.cpp, $(wildcard $(SRC)/*))

#List of project source and header files from each subdirectory
SRCS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
HDRS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.hpp))

#Object file paths to be compiled
OBJS := $(patsubst %.cpp, $(OBJ)/%.o, $(notdir $(SRCS)))

#Compile rule for all o files to output file
$(BIN): $(OBJS)
	$(CC) $(FLAGS) $(OBJS) -o $@


#Compile rule for all cpp files to o files in every subdirectory
.SECONDEXPANSION:

$(OBJ)/%.o: $$(foreach dir,$(DIRS),$$(wildcard $$(dir)/$$*.cpp))
	$(CC) $(FLAGS) -c $^ -o $@

#Clean all o files and binaries
clean:
	rm $(OBJ)/*.o $(BIN) $(RES)/*.csv $(RES)/*.pdb

clean_results:
	rm  $(RES)/*.csv

print_srcs:
	@echo $(SRCS)

print_hdrs:
	@echo $(HDRS)

print_dirs:
	@echo $(DIRS)

print_objs:
	@echo $(OBJS)
