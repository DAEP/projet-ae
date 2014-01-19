#
# Makefile
#

#-------------------------------------------------------------------------------
#		Programs
#-------------------------------------------------------------------------------
CC ?= gcc
CXX ?= g++

#-------------------------------------------------------------------------------
#		Flags
#-------------------------------------------------------------------------------
DEBUG += -g
OPTIMIZATION += #-O3
INCLUDES += 
LIBS += -lm

CFLAGS += -Wall -pedantic -std=c99
CFLAGS += $(DEBUG) $(OPTIMIZATION) $(INCLUDES)
LDFLAGS += $(DEBUG) $(OPTIMIZATION) $(LIBS)

#-------------------------------------------------------------------------------
#		Files
#-------------------------------------------------------------------------------
VPATH = src
OBJ = obj
BIN = bin

C_SRC = main.c problem.c visual.c
C_OBJECTS = $(patsubst %.c,$(OBJ)/%.o,$(C_SRC))

PROGRAM_NAME = code_C_seq
PROGRAM = $(patsubst %,$(BIN)/%,$(PROGRAM_NAME))

#-------------------------------------------------------------------------------
#		Commands
#-------------------------------------------------------------------------------
quiet_cmd_cc_o_c = "  [CC] $@"
      cmd_cc_o_c = $(CC) $(CFLAGS) -c -o $@ $<

quiet_cmd_cc = "  [LD] $@"
      cmd_cc = $(CC) -o $@ $^ $(LDFLAGS)

#-------------------------------------------------------------------------------
#		Rules
#-------------------------------------------------------------------------------
all: $(PROGRAM)

$(C_OBJECTS): | $(OBJ)
$(PROGRAM): | $(BIN)

$(OBJ):
	@mkdir -p $@

$(BIN):
	@mkdir -p $@

$(PROGRAM): $(C_OBJECTS)
	@echo $(quiet_cmd_cc)
	@$(cmd_cc)

$(OBJ)/%.o : %.c
	@echo $(quiet_cmd_cc_o_c)
	@$(cmd_cc_o_c)

.PHONY: clean mrproper rebuild valgrind

clean:
	@rm -rf $(OBJ)

mrproper: clean
	@rm -rf $(PROGRAM)

rebuild: mrproper all

valgrind:
	valgrind -v --leak-check=full --track-origins=yes --show-reachable=yes ./$(PROGRAM)
