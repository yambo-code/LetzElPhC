#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# 1. Locate clang-tidy
# Default to the Homebrew path you provided if it exists, otherwise rely on PATH.
CLANG_TIDY_PATH="/opt/homebrew/Cellar/llvm/21.1.8/bin/clang-tidy"
MPICC=mpicc
#
#
#
#
if [ -x "$CLANG_TIDY_PATH" ]; then
    CLANG_TIDY="$CLANG_TIDY_PATH"
elif command -v clang-tidy &> /dev/null; then
    CLANG_TIDY="clang-tidy"
else
    echo "Error: clang-tidy not found. Please install llvm or add it to your PATH."
    exit 1
fi

# 2. Get MPI include paths
# This queries the mpicc compiler wrapper to find where mpi.h is located.
if command -v $MPICC &> /dev/null; then
    MPI_FLAGS=$($MPICC -show 2>/dev/null | grep -o '\-I[^ ]*' | tr '\n' ' ')
else
    echo "Warning: mpicc not found in PATH. MPI headers might not be resolved."
    MPI_FLAGS=""
fi

# 3. Get Project Specific Includes and Defines from src/make.inc
MAKE_INC="src/make.inc"
EXTRA_INCLUDES=""
if [ -f "$MAKE_INC" ]; then
    # Strip comments (#...), then extract any -I or -D flags (e.g., -I/opt/..., -DCOMPILE_ELPH_DOUBLE)
    EXTRA_INCLUDES=$(sed 's/#.*//' "$MAKE_INC" | grep -o '\-[ID][^ "]*' | tr '\n' ' ')
else
    echo "Warning: $MAKE_INC not found."
fi

# 4. Gather files to analyze
# If you pass a specific file (e.g., ./run-clang-tidy.sh src/main_elphC.c), it analyzes just that.
# Otherwise, it finds all .c files in the src/ directory.
if [ "$#" -gt 0 ]; then
    FILES="$@"
else
    FILES=$(find src -type f -name "*.c")
fi

echo "Using clang-tidy   : $CLANG_TIDY"
echo "Compiler includes  : -Isrc -I. $MPI_FLAGS $EXTRA_INCLUDES"
echo "------------------------------------------------------------------------"

# 5. Run clang-tidy
for file in $FILES; do
    echo "Analyzing $file..."
    # We use '|| true' so the script continues to the next file even if warnings are found
    $CLANG_TIDY "$file" -- -Isrc -I. $MPI_FLAGS $EXTRA_INCLUDES || true
done

echo "------------------------------------------------------------------------"
echo "Analysis complete!"

