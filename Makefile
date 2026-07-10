BUILD_DIR := build
CMAKE ?= cmake

ifeq ($(OS),Windows_NT)
CMAKE_GENERATOR := MinGW Makefiles
CONFIGURE_ARGS := -G "$(CMAKE_GENERATOR)" -DCMAKE_C_COMPILER=gcc
APP := $(BUILD_DIR)/mag_callib_static.exe
else
CMAKE_GENERATOR := Unix Makefiles
CONFIGURE_ARGS := -G "$(CMAKE_GENERATOR)"
APP := $(BUILD_DIR)/mag_callib_static
endif

.PHONY: all configure build clean run

all: build run

configure:
	$(CMAKE) -S . -B $(BUILD_DIR) $(CONFIGURE_ARGS)

build: configure
	$(CMAKE) --build $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

run:
	$(APP)
