COMPILER      := g++
STANDARD      := -std=c++11

WARNINGS      := -Wall -Wextra -Werror -Wpedantic
OPENMP        := -fopenmp
OPTIMIZATIONS := -O3 -march=native -ffast-math

SOURCE_DIR    := ./source
BUILD_DIR     := ./build
OUTPUT_DIR    := .
TARGET        := $(OUTPUT_DIR)/ebeDREENA

CXXFLAGS      := $(STANDARD) $(WARNINGS) $(OPENMP) $(OPTIMIZATIONS)

SRCS          := $(wildcard $(SOURCE_DIR)/*.cpp)
OBJS          := $(patsubst $(SOURCE_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(OUTPUT_DIR)
	$(COMPILER) $(CXXFLAGS) $(OBJS) -o $(TARGET)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(BUILD_DIR)
	$(COMPILER) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET)