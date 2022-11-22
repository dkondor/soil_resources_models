# demography_pars.r -- process results from the parameter exploration
# of population model variants, create plots

library(ggplot2)
library(metR)

colors = c("#5aa732","#005c96","#e75216", "#009ace","#ffd500",
           "#8671a7", "#e78502","#db6aa2","#007939","#9a2846","#815234","#000000")

setwd('/home/dkondor/CSH/HoloSim/demography_paper/demography_models')


# Volterra-model, logistic resource replenishment
tmp1 = read.csv('pars_summary_bt34.out')

p1 = ggplot(tmp1, aes(x=c/r, y=a/r)) + geom_raster(aes(fill=step1))
p1 = p1 + scale_fill_continuous(type = "viridis", limits = c(0, 0.62))
p1 = p1 + metR::geom_contour2(aes(z=step1, label=stat(level)), skip=0,
                              color="white", breaks=c(0.01, 0.05, 0.1, 0.2, 0.4),
                              size = 0.25, label_size = 2,
                              label.placer = label_placer_random())
p1 = p1 + scale_x_continuous(expand = expansion(0,0), limits = c(0,3.0015))  +
  scale_y_continuous(expand = expansion(0,0), limits = c(0,3.0015))
p1 = p1 + theme_bw(8) + labs(fill="Step")

ggsave('volterra_step.pdf', p1, width=3.6, height=2.2)
ggsave('volterra_step.png', p1, width=3.6, height=2.2, dpi=300)


# Volterra-model, linear resource replenishment
tmp2 = read.csv('pars_summary_bte34.out')

p1 = ggplot(tmp2, aes(x=c/r, y=a/r)) + geom_raster(aes(fill=step1))
p1 = p1 + scale_fill_continuous(type = "viridis", limits = c(0, 0.62))
p1 = p1 + metR::geom_contour2(aes(z=step1, label=stat(level)), skip=0,
            color="white", breaks=c(0.01, 0.05, 0.1, 0.2, 0.3),
            size = 0.25, label_size = 2)
p1 = p1 + scale_x_continuous(expand = expansion(0,0), limits = c(0,1.0005))  +
  scale_y_continuous(expand = expansion(0,0), limits = c(0,1.0005))
p1 = p1 + theme_bw(8) + labs(fill="Step")

ggsave('volterra_lin_step.pdf', p1, width=3.6, height=2.2)
ggsave('volterra_lin_step.png', p1, width=3.6, height=2.2, dpi=300)


# other results, re-do the above figures after reading these files
tmp1 = read.csv('pars_summary_bts23.out')
tmp1 = read.csv('pars_summary_bts22.out')
tmp1 = read.csv('pars_summary_bts2.out')
tmp1 = read.csv('pars_summary_bts24.out')


tmp1 = read.csv('pars_summary_bt33.out')
tmp1 = read.csv('pars_summary_bt32.out')
tmp1 = read.csv('pars_summary_bt3.out')

tmp1 = read.csv('pars_summary_bte33.out')
tmp1 = read.csv('pars_summary_bte32.out')
tmp1 = read.csv('pars_summary_bte3.out')

p1 = ggplot(tmp1) + geom_raster(aes(x=c/r, y=a/r, fill=step1))
p1 = p1 + scale_fill_continuous(type = "viridis")
p1 = p1 + geom_contour(aes(x=c/r, y=a/r, z=step1, color=after_stat(level)))
p1 = p1 + scale_color_continuous(type="viridis")
p1 = p1 + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1))


tmp1$step1l = log10(tmp1$step1)


p1 = ggplot(tmp1, aes(x=c/r, y=a/r)) + geom_raster(aes(fill=step1))
p1 = p1 + scale_fill_continuous(type = "viridis", limits=c(0,NA))
# p1 = p1 + geom_contour(aes(x=c/r, y=a/r, z=step1), color="red",
#                        breaks=c(0.01, 0.05, 0.1, 0.2, 0.4))
p1 = p1 + metR::geom_contour2(aes(z=step1, label=stat(level)), skip=0,
                      color="white", breaks=c(0.01, 0.05, 0.1, 0.2, 0.4))
p1 = p1 + theme_bw(6) + labs(fill="Step")


# example contour for the Volterra model case, chi = 10
z1 = 4*pi*pi/(log(10)**2) + 1
tmp12 = data.frame(ar = seq(0, 3, 0.01))
tmp12$cr = (1 + sqrt(1 + z1*tmp12$ar))/z1

p1 = p1 + geom_line(data=tmp12, aes(x=cr, y=ar), color="black")
p1 = p1 + geom_contour(aes(x=c/r, y=a/r, z=step1), breaks=c(0.1), color="red")
# note: there is a mismatch, but not very large
# this can be expected since the estimate based on the eigenvalues
# is not exact


tmp1 = read.csv('pars_summary_bt23.out')
tmp1 = read.csv('pars_summary_bt22.out')
tmp1 = read.csv('pars_summary_bt2.out')

p1 = ggplot(tmp1) + geom_raster(aes(x=c/d, y=b/d-1, fill=step1))
p1 = ggplot(tmp1) + geom_raster(aes(x=c/d, y=b/d-1, fill=overshoot))
p1 = p1 + scale_fill_continuous(trans="log10", limits=c(1e-2,NA))


# plot for the soil-farmer model, new results, increased range of
# parameters
tmp1$step1[is.nan(tmp1$step1)] = 0
p1 = ggplot(tmp1, aes(x=c/r, y=a/r)) + geom_raster(aes(fill=step1))
p1 = p1 + scale_fill_continuous(type = "viridis", option = "A", limits=c(0,NA))
p1 = p1 + theme_bw(8) + labs(fill="Step")
p1 = p1 + scale_x_continuous(expand = expansion(0,0), breaks = c(0,1,2),
                             limits = c(-0.005,2.005))  +
  scale_y_continuous(expand = expansion(0,0), limits = c(-0.05,20.05))
ggsave("soil_step.pdf",p1,width=3.6,height=2.2)
ggsave("soil_step.png",p1,width=3.6,height=2.2,dpi=200)


