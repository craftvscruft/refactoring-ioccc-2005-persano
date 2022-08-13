#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

function runCase() {
  P=$1
  Q=$2
  ./persano $P $Q > out/$P-$Q.gif
  echo Comparing out/$P-$Q.gif
  diff out/$P-$Q.gif golden/$P-$Q.gif
}

runCase 1 1
runCase 3 2
runCase 5 2
runCase 7 2
runCase 4 3
runCase 3 16