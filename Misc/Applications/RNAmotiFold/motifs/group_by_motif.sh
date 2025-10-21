#!/bin/bash

# Usage:
#   ./group_by_motif.sh input.csv group_name output_file.csv


if [ $# -ne 3 ]; then
  echo "Usage: $0 input.csv group_name output_file.csv"
  echo "Example: $0 data.csv fruit grouped/fruit_data.csv"
  exit 1
fi

INPUT_FILE="$1"
GROUP_NAME="$2"
OUTPUT_PATH="$3"

# Check input file
if [ ! -f "$INPUT_FILE" ]; then
  echo "Error: Input file '$INPUT_FILE' not found."
  exit 1
fi

# Extract directory path from output path
OUTPUT_DIR=$(dirname "$OUTPUT_PATH")

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
  echo "Output directory '$OUTPUT_DIR' does not exist. Creating it..."
  mkdir -p "$OUTPUT_DIR"
  if [ $? -ne 0 ]; then
    echo "Failed to create directory: $OUTPUT_DIR"
    exit 1
  fi
fi

# Clear the output file (overwrite if exists)
: > "$OUTPUT_PATH"

# Extract lines matching the group name
awk -F',' -v group="$GROUP_NAME" '
{
    gsub(/^ +| +$/, "", $1)
    gsub(/^ +| +$/, "", $2)
    if ($2 == group) {
        print $1 "," $2
    }
}
' "$INPUT_FILE" > "$OUTPUT_PATH"

echo "Group \"$GROUP_NAME\" written to: $OUTPUT_PATH"
