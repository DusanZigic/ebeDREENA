#!/bin/sh

set -e

COMPILER="g++"
STANDARD="-std=c++11"
WARNINGS="-Wall -Wextra -Werror -Wpedantic"
OPENMP="-fopenmp"
OPTIMIZATIONS="-O3 -march=native -ffast-math"
SOURCE_DIR="./source/"
OUTPUT_DIR="."
TARGET="$OUTPUT_DIR/ebeDREENA"

mkdir -p "$OUTPUT_DIR"

$COMPILER $STANDARD $WARNINGS $OPENMP $OPTIMIZATIONS $SOURCE_DIR*.cpp -o "$TARGET"