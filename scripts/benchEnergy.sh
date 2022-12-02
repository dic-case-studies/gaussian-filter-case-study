#!/usr/bin/env bash

set -e
set -x


declare -a METHODS=(GOLDEN SSE AVX)

host=$1
image=$2

mkdir -p stat/$host/energy
rm -f stat/$host/energy/$image-stats.csv

for method in ${METHODS[@]}
do
    rm -f build/main stat/$host/energy/$image-$method-result.txt

    make build/main OPT="-DNGOLDEN -DNSSE -DNAVX -UASSERT -UN$method"

    echo "Running main ${method} " >> stat/$host/energy/$image-$method-result.txt 
    perf stat -e power/energy-pkg/ ./build/main $image &>> stat/$host/energy/$image-$method-result.txt 
    echo "----------------------------" >> stat/$host/energy/$image-$method-result.txt 
done


AWK=awk
function pattern() {
  id=$1
  label=$2
  PATTERN="
  /energy-pkg/ {
    energy = \$1;
  }
  /seconds time elapsed/ {
    time = \$1;
    printf(\"%s %s %s %s\n\", $id, \"$label\", energy, time);
  }
  "
  echo $PATTERN
}

for ((i=0 ; i < 3; i++))
do
    method=${METHODS[i]}
    PATTERN=`pattern $i ${METHODS[i]}` 
    cat stat/$host/energy/$image-$method-result.txt | $AWK "$PATTERN" >> stat/$host/energy/${image}-stats.csv
done

  echo "                                            \
    reset;                                          \
    set terminal png enhanced large; \
                                                          \
    set title \"Gaussian Filter Energy Benchmark\";                        \
    set xlabel \"Methods\";                             \
    set ylabel \"Joules\";                     \
    set yrange [0:150];                 \
    unset key;                                      \
    set boxwidth 0.5;                                       \
    set style fill solid;                                \
                                                          \
    plot \"stat/$host/energy/${image}-stats.csv\" using 1:3:xtic(2) with boxes
" | gnuplot > stat/$host/energy/${image}-performance.png
