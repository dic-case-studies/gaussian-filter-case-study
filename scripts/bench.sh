#!/usr/bin/env bash

set -e
set -x

host=$1
image=$2
mkdir -p stat/$host

rm -rf stat/$host/${image}-result.txt stat/$host/${image}-stats.csv out/*

./build/main $image > stat/$host/${image}-result.txt

METHOD=AVX

processor=$(uname -m)

if [ $processor = "arm64" ]
then
    METHOD=NEON 
fi

declare -a methods=(GOLDEN SSE $METHOD)


AWK=awk
function pattern() {
  id=$1
  label=$2
  PATTERN="
  /Time elapsed $label/ {
    time = \$(NF-1);
    printf(\"%s %s %s\n\", $id, \"$label\", time);
  }
  "
  echo $PATTERN
}

for ((i=0 ; i < 3; i++))
do
    PATTERN=`pattern $i ${methods[i]}` 
    cat stat/$host/${image}-result.txt | $AWK "$PATTERN" >> stat/$host/${image}-stats.csv
done

echo "                                            \
reset;                                          \
set terminal png enhanced large; \
                                                        \
set title \"Gaussian Filter Benchmark\";                        \
set xlabel \"METHODS\";                             \
set ylabel \"Execution time(s)\";                     \
set yrange [0:15];                 \
unset key;                                      \
set boxwidth 0.5;                                       \
set style fill solid;                                \
                                                        \
plot \"stat/$host/${image}-stats.csv\" using 1:3:xtic(2) with boxes; \
" | gnuplot > stat/$host/${image}-performance.png