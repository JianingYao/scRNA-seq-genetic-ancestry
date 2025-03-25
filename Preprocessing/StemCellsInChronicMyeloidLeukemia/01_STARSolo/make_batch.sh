#!/bin/bash

source_dir="/scratch1/yaojiani/StemCellsInChronicMyeloidLeukemia"
base_dir="/scratch1/yaojiani/StemCellsInChronicMyeloidLeukemia/batches"

for i in {1..2}; do
  mkdir -p "$base_dir/batch$i"
done

batch=1

for file in "$source_dir"/*; do
  mv "$file" "$base_dir/batch$batch"
  batch=$((batch % 2 + 1))
done
