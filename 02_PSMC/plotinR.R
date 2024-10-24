library(ggplot2)
library(tidyverse)

## FINAL

## Read in the data:
dat <- read.csv ("PSMC_Dec27_done_noQinghai.csv", header=TRUE)
## Grab only past 10K years ago
dat2 <- subset(dat, Time > 9000)

############################################################################
### Plotting PSMC with other gray wolf individuals to check parameters #####
############################################################################
## Plot correction
p <- ggplot(data=dat2, aes(x=Time, y=Ne, color=Population)) +
  geom_step(size=1, alpha=0.7) + theme_classic() + scale_x_log10() + scale_color_manual(values=c("maroon4", "darkseagreen4", "orange2", "darkseagreen4", "yellowgreen", "blue", "yellowgreen", "darkseagreen4", "green2", "maroon4", "lightseagreen", "lightseagreen", "lightseagreen", "lightseagreen", "lightseagreen", "lightseagreen", "lightseagreen", "lightseagreen", "yellowgreen", "yellowgreen", "yellowgreen", "yellowgreen", "darkseagreen4", "lightseagreen", "green2", "yellowgreen", "yellowgreen", "yellowgreen", "yellowgreen", "yellowgreen", "yellowgreen", "yellowgreen", "yellowgreen", "orange2", "orange2" ))
p + 
  facet_grid(Sample ~ .)

ggsave("July17_PSMC_testparameters.tiff", width=6,height=5)
