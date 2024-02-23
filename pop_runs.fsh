#!/usr/bin/fish
# pop_runs.fsh -- run population models with noise
# This script uses the Fish shell syntax, see here:
# https://fishshell.com/
# It is assumed that this script is run from the directory where the code is.

# Compile the code before running:
# g++ -o pm1 pop_model1.cpp -lgsl -lgslcblas -O3 -march=native -lm


# 1. Main soil model with two parameter combinations and noise (repeated realizations)

set outdir model_runs
set outdir2 $outdir/rep1
mkdir -p $outdir2

# 1.1. Early agriculturalists case (Fig. 3)
for i in (seq 100)
./pm1 -t LL -a 0.3 -c 0.2731 -g 1 -r 0.02 -K 10000 -x 100 -y 10000 -T 10000 -s $i -N 0.05 > $outdir2/LL1N5_$i.out
./pm1 -t L  -a 0.3 -c 0.2    -g 1 -r 0.02 -K 10000 -x 100 -y 10000 -T 10000 -s $i -N 0.05 >  $outdir2/L1N5_$i.out
end

# 1.2. Slow resource regrowth case (Fig. 4)
for i in (seq 100)
./pm1 -t LL -a 0.1 -c 0.06 -g 1 -r 0.02 -K 10000 -x 100 -y 10000 -T 10000 -N 0.05 -s $i > $outdir2/LL2N5_$i.out
./pm1 -t L  -a 0.1 -c 0.03835 -g 1 -r 0.02 -K 10000 -x 100 -y 10000 -T 10000 -s $i -N 0.05  >  $outdir2/L2N5_$i.out
end



# 2. Turchin-Korotayev warfare model, example (Fig. 5)
./pm1 -t T -a 0.005 -b 0.004 -r 0.02 -c 1 -K 10000 -T 2000 -Nx 0.05 -s 1 -x 1000 -y 0.01 > $outdir/tk1.out



