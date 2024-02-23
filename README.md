# soil_resources_models

Code and scripts for numerical solutions of the models in our paper:<br>
Kondor D, Turchin P (2022). Soil fertility depletion is not a credible mechanism for population boom/bust cycles in agricultural societies. https://osf.io/preprints/socarxiv/rjns2/

This repository includes C++ code to calculate numerical solutions of resource-consumer models. It also includes shell scripts to calculate solutions with the parameters combinations used in our paper and R and Gnuplot scripts that create the figures in the paper.

Requirements:
 - [GSL](https://www.gnu.org/software/gsl/)
 - [Gnuplot](http://gnuplot.info)
 - [R](https://www.r-project.org/) with [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)

Before running the scripts, compile the code:
```
g++ -o pm1 pop_model1.cpp -lgsl -lgslcblas -O3 -march=native -lm
```

