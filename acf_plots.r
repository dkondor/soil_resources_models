library(ggplot2)

# !! TODO: set the correct working directory !!
# setwd('/home/dkondor/CSH/HoloSim/demography_paper/demography_models/soil_resources_models/')

# maximum lag for ACFs
lmax = 1000
# base directory with simulation output
outdir = c('model_runs/')

# analyse repeated runs
acf_repeat = function(fn1, N) {
  acf1 = data.frame()
  for(i in 1:N) {
    tmp1 = read.table(paste0(fn1, i, ".out"), header = FALSE)
    names(tmp1) = c('t', 'pop', 'res')
    tmp1 = tmp1[tmp1$t >= 1000,]
    
    acfp1 = acf(tmp1$pop, lag.max = lmax, plot = FALSE)
    acfp1 = data.frame(lag = acfp1$lag, pop = acfp1$acf)
    acfp1$i = i
    acf1 = rbind(acf1, acfp1)
  }
  return(acf1)
}

# create figures comparing the biotic and abiotic model variants
create_fig = function(fn1, fn2) {
  acf1 = acf_repeat(fn1, 100)
  acf2 = acf_repeat(fn2, 100)
  acf1$Resources = 'Biotic'
  acf2$Resources = 'Abiotic'
  acf3 = rbind(acf1, acf2)
  
  p1 = ggplot(mapping=aes(x=lag, y=pop, group=interaction(Resources,i), color=Resources))
  p1 = p1 + geom_line(data=acf3, size = 0.1)
  p1 = p1 + scale_color_manual(values=c("#e75216", "#005c96"))
  p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("ACF")
  
  # optional: plot the mean of all ACFs
  # acf1m = aggregate(acf1[c("pop", "i")], by=acf1["lag"], FUN=mean)
  # acf2m = aggregate(acf2[c("pop", "i")], by=acf2["lag"], FUN=mean)
  
  # p1 + geom_line(data=acf1m) + geom_line(data=acf2m)
  return(p1)
}

# 1.1. Early agriculturalists case (Fig. 3, right panel)
fn1 = paste0(outdir, "rep1/LL1N5_")
fn2 = paste0(outdir, "rep1/L1N5_")
p1 = create_fig(fn1, fn2)
fn1 = paste0(outdir, 'acf_cmp1')
ggsave(paste0(fn1, '.pdf'), p1, width=3.4, height=2)
ggsave(paste0(fn1, '.png'), p1, width=3.4, height=2, dpi=300)

# 1.2. Slow resource regrowth case (Fig. 4, right panel)
fn1 = paste0(outdir, "rep1/LL2N5_")
fn2 = paste0(outdir, "rep1/L2N5_")
p1 = create_fig(fn1, fn2)
fn1 = paste0(outdir, 'acf_cmp2')
ggsave(paste0(fn1, '.pdf'), p1, width=3.4, height=2)
ggsave(paste0(fn1, '.png'), p1, width=3.4, height=2, dpi=300)

