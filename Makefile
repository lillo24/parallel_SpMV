# Made fully with AI not gonna lie, I didn't do Operating System




# Makefile
# Build:    make
# Run:      make run FILE=1138_bus.mtx
# Clean:    make clean
# Debug:    make DEBUG=1

# --- toolchain & flags ---
CC      ?= gcc
OPT     := -O3
WARN    := -Wall -Wextra
STD     := -std=c11
ARCH    := -march=native
OPENMP  := -fopenmp

# switch to debug with: make DEBUG=1
ifdef DEBUG
  OPT := -Og -g
endif

CFLAGS_COMMON := $(STD) $(OPT) $(ARCH) $(OPENMP)
CFLAGS_MAIN   := $(CFLAGS_COMMON) $(WARN)
CFLAGS_MMIO   := $(CFLAGS_COMMON) -w   # silence 3rd-party mmio.c warnings
LDFLAGS       := $(OPENMP)

# --- targets ---
TARGET  := mmread
OBJS    := main.o mmio.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

# compile your program (keeps warnings on)
main.o: main.c mmio.h
	$(CC) $(CFLAGS_MAIN) -c $<

# compile mmio.c (quiet; it's third-party)
mmio.o: mmio.c mmio.h
	$(CC) $(CFLAGS_MMIO) -c $<

# convenience runner: make run FILE=path/to/file.mtx
FILE ?= 1138_bus.mtx
run: $(TARGET)
	./$(TARGET) $(FILE)

clean:
	rm -f $(TARGET) *.o
.PHONY: all run clean