# pop_plots.gp -- Gnuplot commands to create figures with example behavior of models
outdir = 'model_runs'
outdir2 = outdir.'/rep1'

# figure setup
se xl 'Time [years]'
se yl 'Population [x1000]' off 1
se ytics 2
se xr [0:2000]


# 1.1. Early agriculturalists case (Fig. 3, left panel)
  p outdir2.'/LL1N5_1.out' u 1:($2/1000) w l lw 2 lc rgbcolor '#005c96' not
rep outdir2.'/L1N5_1.out'  u 1:($2/1000) w l lw 2 lc rgbcolor '#e75216' not
se te post eps color solid size 2.8,2
se out outdir.'/res_cmp1.eps'
rep
se out
se te wxt

# 1.2. Slow resource regrowth case (Fig. 4, left panel)
  p outdir2.'/LL2N5_1.out' u 1:($2/1000) w l lw 2 lc rgbcolor '#005c96' not
rep outdir2.'/L2N5_1.out' u 1:($2/1000) w l lw 2 lc rgbcolor '#e75216' not
se te post eps color solid size 2.8,2
se out outdir.'/res_cmp2.eps'
rep
se out
se te wxt


# 2. Turchin-Korotayev warfare model (Fig. 5)
p outdir.'/tk1.out' u 1:($2/1000) w l lw 2 lc -1 not

se te post eps color solid size 2.8,2
se out outdir.'/res_tk1.eps'
rep
se out
se te wxt



