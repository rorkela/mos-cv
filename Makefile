# Compiler and flags
CC      := gcc
CFLAGS  := -Wall -Wextra -Ofast -I./src 
LDFLAGS := -lm -fopenmp

# Target executable
TARGET  := build/moscv
CONFIG := build/config.txt
# Source files (recursive)
SRCS := \
    src/main.c \
    src/carrier_continuity/carrier.c \
    src/fileio/fileio.c \
    src/parameter_fetch/parameter_fetch.c \
    src/poisson/poisson.c \
    src/solve_c/solve_c.c

# Object files
OBJS := $(SRCS:src/%.c=build/%.o)

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

# Compile
build/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)



# Config file generation
$(CONFIG):$(TARGET)
	
config:$(CONFIG)
clean:
	rm -rf build

# Phony targets
.PHONY: all clean
