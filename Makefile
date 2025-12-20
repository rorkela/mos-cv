.PHONY: all configure build clean rebuild

BUILD_DIR := build

all: build

configure:
	@mkdir -p $(BUILD_DIR)
	cmake -S . -B $(BUILD_DIR)

build: configure
	cmake --build $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

rebuild: clean all
