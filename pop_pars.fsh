#!/usr/bin/fish
# pop_pars.fsh -- explore the parameter space of population models


##########################################################################################
# base runs, loops are in the C++ code
# 1. Volterra model
./pm1 -t bt -d 0.01 -Mb 0.01005 0.0001 0.04 -Mc 0.00005 0.0001 0.03 -a - -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bt2.out
./pm1 -t bt -d 0.01 -Mb 0.0101 0.0002 0.11 -Mc 0.0001 0.0002 0.1 -a - -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bt22.out
./pm1 -t bt -d 0.01 -Mb 0.0101 0.0002 0.11 -Mc 0.0001 0.00002 0.01 -a - -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bt23.out

# 2. linear resource growth
./pm1 -t bte -d 0.01 -Mb 0.01005 0.0001 0.04 -Mc 0.00005 0.0001 0.03 -a - -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bte2.out
./pm1 -t bte -d 0.01 -Mb 0.0101 0.0002 0.11 -Mc 0.0001 0.0002 0.1 -a - -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bte22.out
./pm1 -t bte -d 0.01 -Mb 0.0101 0.0002 0.11 -Mc 0.0001 0.00002 0.01 -a - -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bte23.out

# 3. new soil model
./pm1 -t s2 -r 0.01 -Ma 0.00005 0.0001 0.03 -Mc 0.00005 0.0001 0.03 -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bts2.out
./pm1 -t s2 -r 0.01 -Ma 0.0001 0.0002 0.1 -Mc 0.0001 0.0002 0.1 -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bts22.out
./pm1 -t s2 -r 0.01 -Ma 0.0001 0.0002 0.1 -Mc 0.0001 0.00002 0.01 -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bts23.out
# new version, a up to 20, c up to 2
./pm1 -t s2 -r 0.01 -Ma 0.0001 0.0005 0.2 -Mc 0.0001 0.00005 0.02 -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bts24.out


###############################################################################################
# new cases for the Volterra model, use a parametrization more similar to the other case
./pm1 -t bt -r 0.01 -Ma 0.00005 0.0001 0.03 -Mc 0.00005 0.0001 0.03 -x 40 -y 10000 -K 10000 -Cs 2 -CS -T 1000000 > pars_summary_bt3.out
./pm1 -t bt -r 0.01 -Ma 0.0001 0.0002 0.1 -Mc 0.0001 0.0002 0.1 -x 40 -y 10000 -K 10000 -Cs 2 -CS -T 100000 > pars_summary_bt32.out
./pm1 -t bt -r 0.01 -Ma 0.0001 0.0002 0.1 -Mc 0.0001 0.00002 0.01 -x 40 -y 10000 -K 10000 -Cs 2 -CS -T 100000 > pars_summary_bt33.out

./pm1 -t bt -r 0.01 -Ma 0.00003 0.00003 0.03 -Mc 0.00003 0.00003 0.03 -x 40 -y 10000 -K 10000 -Cs 2 -CS -T 1000000 > pars_summary_bt34.out

./pm1 -t bte -r 0.01 -Ma 0.00005 0.0001 0.03 -Mc 0.00005 0.0001 0.03 -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bte3.out
./pm1 -t bte -r 0.01 -Ma 0.0001 0.0002 0.1 -Mc 0.0001 0.0002 0.1 -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bte32.out
./pm1 -t bte -r 0.01 -Ma 0.0001 0.0002 0.1 -Mc 0.0001 0.00002 0.01 -x 40 -y 10000 -K 10000 -Cs 2 -CS > pars_summary_bte33.out

./pm1 -t bte -r 0.01 -Ma 0.00001 0.00001 0.01 -Mc 0.00001 0.00001 0.01 -x 40 -y 10000 -K 10000 -Cs 2 -CS -T 1000000 > pars_summary_bte34.out




