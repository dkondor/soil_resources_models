#!/usr/bin/fish
# pop_runs.fsh -- run population models with and without noise
# This script uses the Fish shell syntax, see here:
# https://fishshell.com/
# It also uses the runmultiple utility from here:
# https://github.com/dkondor/runmultiple
# (assumed to be installed as ~/program/misc_utils/runmultiple/rm)

# output directory
mkdir noise_res

# four versions:
# A version: a = 0.4, c = 0.2, r = 0.02 (Neolithic farmers' case)
# B version: a = 0.03, c = 0.02, r = 0.02 (salinization case, short timescales)
# C version: a = 0.002, c = 0.002, r = 0.02 (salinization case, longer time interval)
# D version: a = 0.06, c = 0.02, r = 0.03 (no specific case, all models have pseudo-oscillations)

set models 1 2 3 4
set models2 bt bte s2 t
set par_ix 1 1 1 2
set min_pars 1 5
set max_pars 4 5
set bpar "-b -" "-b -" "" ""
set parnames A B C D ""
set pars "-a 0.4 -c 0.2 -r 0.02" "-a 0.03 -c 0.02 -r 0.02" "-a 0.002 -c 0.002 -r 0.02" "-a 0.06 -c 0.02 -r 0.03" "-a 0.005 -b 0.004 -r 0.02 -c 1"
set xy0 "-x 40 -y 10000" "-x 40 -y 0"


for m in $models
for p in (seq $min_pars[$par_ix[$m]] $max_pars[$par_ix[$m]])
echo ./pm1 -t $models2[$m] $pars[$p] $bpar[$m] $xy0[$par_ix[$m]] -K 10000 -T 10000 \> noise_res/res0_$models[$m]$parnames[$p].out
end
end | ~/program/misc_utils/runmultiple/rm -t 1


for NN in x y
for n in 0.01 0.02 0.03 0.04 0.05 0.1
for m in $models
for p in (seq $min_pars[$par_ix[$m]] $max_pars[$par_ix[$m]])
echo ./pm1 -t $models2[$m] $pars[$p] $bpar[$m] $xy0[$par_ix[$m]] -K 10000 -T 20000 -N$NN $n -s (head -c 4 /dev/urandom | hexdump -e "\"%u\"") \> noise_res/noise_$models[$m]$parnames[$p]_$NN$n.out
end
end
end
end | ~/program/misc_utils/runmultiple/rm -t 4






