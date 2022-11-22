# noise_plots.r -- plots of model behavior with noise for the main text and SI

library(ggplot2)

# note: set to the correct base directory
setwd('/home/dkondor/CSH/HoloSim/demography_paper/demography_models')

base_dir = 'noise_res/'
roman14_1 = function(x) {
  return(switch(x, 'I', 'II', 'III', 'IV')) 
}
roman14 = function(x) {
  return(sapply(x, roman14_1))
}


base_fn = 'noise_'
fn1 = c('A', 'B', 'C', 'D')

noise_levels = c(0.01, 0.02, 0.03, 0.04, 0.05) # extended version
var = 'x' # original case

all_res = data.frame()

for(i in 1:4) {
  j2 = if(i < 4) length(fn1) else 1
  for(j in 1:j2) {
    for(n in noise_levels) {
      fn2 = if(i < 4) fn1[j] else ''
      tmp1 = read.table(paste0(base_dir, base_fn, i, fn2, '_', var, n, '.out'),
                        header = FALSE)
      names(tmp1) = c('t', 'pop', 'res')
      tmp1$model = i
      tmp1$var = fn2
      tmp1$n = n
      all_res = rbind(all_res, tmp1)
    }
  }
}


h1 = 6
h2 = 3
suffix = '_1'
noise_levels2 = noise_levels

# plot of all variants and all noise levels
for(x in fn1) {
  tmp1 = all_res[all_res$var == x,]
  tmp1 = tmp1[tmp1$n %in% noise_levels2,]
  tmp1 = tmp1[tmp1$t >= 10000,]
  tmp1$t = tmp1$t - 10000
  tmp1$model2 = paste0("Model ", roman14(tmp1$model))
  tmp1$n2 = paste0("Noise level: ", 100*tmp1$n, "%")
  p1 = ggplot(tmp1) + geom_line(aes(x=t, y=pop), size=0.3)
  p1 = p1 + facet_grid(n2~model2)
  p1 = p1 + scale_x_continuous(limits=c(0,6000), breaks = c(0,2000,4000,6000))
  p1 = p1 + scale_y_continuous(limits=c(0,NA))
  p1 = p1 + theme_bw(8) + xlab("Time [years]") + ylab("Population")
  ggsave(paste0(base_dir, 'model', x, suffix, '.pdf'), p1, width=6.8, height=h1)
  ggsave(paste0(base_dir, 'model', x, suffix, '.png'), p1, width=6.8, height=h1, dpi=200)
}

# separately for model IV
tmp1 = all_res[all_res$model == 4,]
tmp1 = tmp1[tmp1$n %in% noise_levels2,]
tmp1 = tmp1[tmp1$t >= 10000,]
tmp1$t = tmp1$t - 10000
tmp1$n2 = paste0("Noise level: ", 100*tmp1$n, "%")
p1 = ggplot(tmp1) + geom_line(aes(x=t, y=pop), size=0.3)
p1 = p1 + facet_wrap(~n2)
p1 = p1 + scale_x_continuous(limits=c(0,6000), breaks = c(0,2000,4000,6000))
p1 = p1 + scale_y_continuous(limits=c(0,NA))
p1 = p1 + theme_bw(8) + xlab("Time [years]") + ylab("Population")
ggsave(paste0(base_dir, 'model', 4, suffix, '.pdf'), p1, width=6.8, height=h2)
ggsave(paste0(base_dir, 'model', 4, suffix, '.png'), p1, width=6.8, height=h2, dpi=200)


# re-do the above, only using noise levels of 1%, 3%, 5%
noise_levels2 = c(0.01, 0.03, 0.05)
suffix = "_2"
h1 = 4.5
h2 = 1.8


# figure for the main paper -- only 5% noise, three model and three
# parameter combinations
tmp1 = all_res[all_res$n == 0.05,]
tmp1 = tmp1[tmp1$model != 4,]
tmp1 = tmp1[tmp1$var %in% c('A', 'B', 'C'),]
tmp1 = tmp1[tmp1$t >= 10000,]
tmp1$t = tmp1$t - 10000
tmp1$model2 = paste0("Model ", roman14(tmp1$model))
eqnum = function(x) {
  return(switch(x, "Volterra model (logistic)", "Volterra model (abiotic)", "Agriculture model"))
}
tmp1$eq = sapply(tmp1$model, eqnum)
tmp1$eq = factor(tmp1$eq, levels = c("Agriculture model", "Volterra model (logistic)", "Volterra model (abiotic)"))
tmp1$var2 = paste0("Case ", tmp1$var)

p1 = ggplot(tmp1) + geom_line(aes(x=t, y=pop), size=0.3)
p1 = p1 + facet_grid(var2~eq)
p1 = p1 + scale_x_continuous(limits=c(0,6000), breaks = c(0,2000,4000,6000))
p1 = p1 + scale_y_continuous(limits=c(0,20700))
p1 = p1 + theme_bw(8) + xlab("Time [years]") + ylab("Population")
ggsave(paste0(base_dir, 'noise_cmp', '.pdf'), p1, width=6.8, height=4.5)
ggsave(paste0(base_dir, 'noise_cmp', '.png'), p1, width=6.8, height=4.5, dpi=200)

# only cases B and C (Fig. 4)
tmp1 = tmp1[tmp1$var %in% c('B', 'C'),]
ggsave(paste0(base_dir, 'noise_cmp2BC', '.pdf'), p1, width=6.8, height=4)
ggsave(paste0(base_dir, 'noise_cmp2BC', '.png'), p1, width=6.8, height=4, dpi=200)


# only case A
tmp1 = tmp1[tmp1$var == 'A',]
p1 = ggplot(tmp1) + geom_line(aes(x=t, y=pop), size=0.3)
p1 = p1 + facet_grid(eq~.)
p1 = p1 + scale_x_continuous(limits=c(0,6000), breaks = c(0,2000,4000,6000))
p1 = p1 + scale_y_continuous(limits=c(0,15000))
p1 = p1 + theme_bw(8) + xlab("Time [years]") + ylab("Population")
ggsave(paste0(base_dir, 'noise_cmp', '.pdf'), p1, width=4.5, height=4.5)
ggsave(paste0(base_dir, 'noise_cmp', '.png'), p1, width=4.5, height=4.5, dpi=200)


# only model 4, only 5% noise
tmp1 = all_res[all_res$n == 0.05,]
tmp1 = tmp1[tmp1$model == 4,]
tmp1 = tmp1[tmp1$t >= 10000,]
tmp1$t = tmp1$t - 10000

p1 = ggplot(tmp1) + geom_line(aes(x=t, y=pop), size=0.3)
p1 = p1 + scale_x_continuous(limits=c(0,6000), breaks = c(0,2000,4000,6000))
p1 = p1 + theme_bw(8) + xlab("Time [years]") + ylab("Population")
ggsave(paste0(base_dir, 'model4_5', '.pdf'), p1, width=3.2, height=2)
ggsave(paste0(base_dir, 'model4_5', '.png'), p1, width=3.2, height=2, dpi=200)



####################################################################
# ACFs with noise
lmax = 3000
acf_all = data.frame()
for(i in 1:4) {
  j2 = if(i < 4) length(fn1) else 1
  for(j in 1:j2) {
    fn2 = if(i < 4) fn1[j] else ''
    tmp1 = all_res[all_res$model == i & all_res$var == fn2,]
    for(n in noise_levels) {
      tmp2 = tmp1[tmp1$n == n & tmp1$t >= 10000,]
      tmp2 = tmp2[order(tmp2$t),]
      acf1 = acf(tmp2$pop, lag.max = lmax, plot = FALSE)
      acf2 = data.frame(lag = acf1$lag, acf = acf1$acf)
      acf2$model = i
      acf2$var = fn2
      acf2$n = n
      acf_all = rbind(acf_all, acf2)
    }
  }
}
 


h1 = 6
h2 = 3
suffix = '_1'
noise_levels2 = noise_levels

# plot of all variants and all noise levels
for(x in fn1) {
  tmp1 = acf_all[acf_all$var == x,]
  tmp1 = tmp1[tmp1$n %in% noise_levels2,]
  tmp1$model2 = paste0("Model ", roman14(tmp1$model))
  tmp1$n2 = paste0("Noise level: ", 100*tmp1$n, "%")
  p1 = ggplot(tmp1) + geom_line(aes(x=lag, y=acf), size=0.3)
  p1 = p1 + facet_grid(n2~model2)
  # add a fake x-axis
  p1 = p1 + geom_line(data=data.frame(x=c(-50,3050),y=c(0,0)),aes(x=x,y=y),size=0.1)
  p1 = p1 + scale_x_continuous(limits=c(-50,3050), expand=c(0,0), breaks = c(0,1000,2000,3000))
  p1 = p1 + theme_bw(8) + xlab("Lag [years]") + ylab("ACF")
  ggsave(paste0(base_dir, 'acf_', x, suffix, '.pdf'), p1, width=6.8, height=h1)
  ggsave(paste0(base_dir, 'acf_', x, suffix, '.png'), p1, width=6.8, height=h1, dpi=200)
}

# separately for model IV
tmp1 = acf_all[acf_all$model == 4,]
tmp1 = tmp1[tmp1$n %in% noise_levels2,]
tmp1$n2 = paste0("Noise level: ", 100*tmp1$n, "%")
p1 = ggplot(tmp1) + geom_line(aes(x=lag, y=acf), size=0.3)
p1 = p1 + facet_wrap(~n2)
# add a fake x-axis
p1 = p1 + geom_line(data=data.frame(x=c(-50,3050),y=c(0,0)),aes(x=x,y=y),size=0.1)
p1 = p1 + scale_x_continuous(limits=c(-50,3050), expand=c(0,0), breaks = c(0,1000,3000))
p1 = p1 + theme_bw(8) + xlab("Lag [years]") + ylab("ACF")
ggsave(paste0(base_dir, 'acf_', 4, suffix, '.pdf'), p1, width=6.8, height=h2)
ggsave(paste0(base_dir, 'acf_', 4, suffix, '.png'), p1, width=6.8, height=h2, dpi=200)


# re-do the above, only using noise levels of 1%, 3%, 5%
noise_levels2 = c(0.01, 0.03, 0.05)
suffix = "_2"
h1 = 4.5
h2 = 1.8


     
##################################################################
# plot of model behavior without noise (SI figs.)
base_fn = 'res0_'
fn1 = c('A', 'B', 'C', 'D')

all_res0 = data.frame()
for(i in 1:4) {
  j2 = if(i < 4) length(fn1) else 1
  for(j in 1:j2) {
    fn2 = if(i < 4) fn1[j] else ''
    tmp1 = read.table(paste0(base_dir, base_fn, i, fn2, '.out'),
                        header = FALSE)
    names(tmp1) = c('t', 'pop', 'res')
    tmp1$model = i
    tmp1$var = fn2
    all_res0 = rbind(all_res0, tmp1)
  }
}

tmp2 = all_res0[all_res0$model != 4,]
tmp2 = tmp2[tmp2$var != 'D',]
tmp2$model2 = paste0("Model ", roman14(tmp2$model))


p1 = ggplot(tmp2) + geom_line(aes(x=t, y=pop, color="Population"))
p1 = p1 + geom_line(aes(x=t, y=res, color="Resources"))
p1 = p1 + facet_grid(var~model2, scales = "free_y")
p1 = p1 + scale_x_continuous(limits=c(0,4000), breaks = c(0,2000,4000))
p1 = p1 + theme_bw(8) + xlab("Time [years]") + ylab("Population")
p1 = p1 + scale_color_manual(values=c("black", "red"))
p1 = p1 + labs(color="")
p1 = p1 + theme(legend.position="bottom")
ggsave(paste0(base_dir, 'res0', '.pdf'), p1, width=6.8, height=5)
ggsave(paste0(base_dir, 'res0', '.png'), p1, width=6.8, height=5, dpi=200)


