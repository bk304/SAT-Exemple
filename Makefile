SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = $(BIN_DIR)/obj
MINI_DIR = minisat

SOLV_EXE = $(BIN_DIR)/solve
JOGO_EXE = $(BIN_DIR)/jogo

SOURCES =                        \
	$(MINI_DIR)/core/Solver.cc   \
	$(MINI_DIR)/utils/Options.cc \
	$(MINI_DIR)/utils/System.cc  \

SOLV_SOURCES = $(SOURCES) $(SRC_DIR)/main.cc
JOGO_SOURCES = $(SOURCES) $(SRC_DIR)/jogo.cc

SOLV_OBJECTS = $(patsubst %.cc,$(OBJ_DIR)/%.o,$(SOLV_SOURCES))
JOGO_OBJECTS = $(patsubst %.cc,$(OBJ_DIR)/%.o,$(JOGO_SOURCES))

$(info $(OBJECTS))

CC = g++
CFLAGS = -I. -I$(MINI_DIR) -I$(SRC_DIR) -fpermissive

all: $(SOLV_EXE) $(JOGO_EXE)

$(SOLV_EXE): $(SOLV_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

$(JOGO_EXE): $(JOGO_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: %.cc
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: all clean
.SUFFIXES: