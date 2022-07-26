#!/bin/bash

set -eu

PREFIX="/home/tom/src/spin-ed/example"

for numLocales in 4 3 2 1; do
  for yamlPath in \
    "heisenberg_chain_10.yaml" \
    "heisenberg_kagome_12.yaml" \
    "heisenberg_kagome_16.yaml" \
    "heisenberg_square_4x4.yaml" \
    "heisenberg_square_5x5.yaml" \
    "heisenberg_triangular_19.yaml"; do

    ./test_basis_construction -nl $numLocales \
      --kInputBasisPath "$PREFIX/$yamlPath" \
      --kInputDataPath "$PREFIX/data/${yamlPath%.yaml}.h5" \
      --kOutputDataPath "output.h5"
    # h5diff "output.h5" "$referenceOutput" \
    #   "/representatives" "/basis/representatives"
  done
done
