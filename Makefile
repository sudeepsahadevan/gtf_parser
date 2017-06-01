# gcc stuff
CC=g++
CFLAGS=-c -Wall -g3 -std=c++11 -O3 
LIBS=-lboost_iostreams
RM=rm -rf
# other stuff
BIN=gtfparser
BINF=bin
SRC=$(wildcard src/*.cpp)
OBJ=$(patsubst src/%,%,$(SRC:%.cpp=%.o))

all: bin	

bin: $(BIN)
	@mkdir -p $(BINF)
	@mv $(BIN) $(BINF)/
	@echo 'executable file: ' $(BINF)/$(BIN)
	@echo 'done...'

$(BIN): $(OBJ)
	@$(CC) $< $(LIBS) -o $@
	@echo 'generating executable'

$(OBJ): $(SRC)
	@echo 'compiling file: ' $<	
	@$(CC) $(CFLAGS) $< -o $@
	
clean:	
	@$(RM) *.o $(BINF)
	@echo 'done cleaning up'
