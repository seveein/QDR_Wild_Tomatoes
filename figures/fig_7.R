# fig. 7
# Modelling AUDPC from lme-coefficients
# 2024_02_01

library(tidyverse)
library(lsmeans)
library(nlme)
rm(list=ls())
setwd("C:/Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/")
data <- read.csv("Spen_subset_Slopes_AUDPC.csv",
                 header=T)

plot(data$audpc)

data$genotype <- as.factor(data$genotype)
data$species <- as.factor(data$species)
data$start_date <- as.factor(data$start_date)

boxplot(audpc ~ genotype, data=data, main="AUDPC", xlab="", ylab="", las=3, 
        col=c("darkgreen", "tomato", "gold"))

# model

mod_contr <- lme(audpc ~  genotype*lag*LDT, data=data,
                 random=~1|start_date,
                 weights = varIdent(form=~1|genotype),
                 control=list(maxIter =500, msMaxIter= 500, 
                              niterEM=500, msMaxEval=500, opt="nlminb"))
mod_contr2 <- update(mod_contr, . ~ 0 + genotype + genotype:lag + genotype:LDT + genotype:lag:LDT)
summary(mod_contr2)

(anova_coefficients<-data.frame(anova(mod_contr)) %>% 
    mutate(p.value = format.pval(pv = p.value,
                                 digits = 2, 
                                 eps = 0.001,
                                 nsmall = 3))
)

(output_lme <- mod_contr2$coefficients$fixed )

# Define ranges for LDT
# hours to LDT(mins) log(x * 60)  
# (8640-1440)*0.01
ldt_range <- seq(5, 7.31, 0.0231)

# define lag-range 
lag_range <- seq(1440, 4320, 28.8)
# define genotypes
genotypes <- c("LA1282", "LA1809", "LA1941")
# Generate empty dataframes to store modeled data 
audpc_values <- data.frame()

# Iterate over LDT and lag values for each genotype

for (genotype in genotypes) {
  for (LDT in ldt_range) {
    for (lag in lag_range) {
      # Calculate AUDPC for each LDT and lag combination
      coefs <- output_lme[str_detect(names(output_lme), genotype)]
      audpc_estimates <- coefs[1] + coefs[2]*lag + coefs[3]*LDT + coefs[4]*lag*LDT
      
      audpc_value_temp <- data.frame(
        genotype=genotype, 
        AUDPC = audpc_estimates, 
        lag=lag, 
        LDT=LDT
      )
      
      # Construct and print the corresponding DataFrame
      audpc_values = rbind(audpc_values,audpc_value_temp)
    }
  }
}


# Viz

(ggplot(audpc_values, aes(x=lag, y=LDT, fill=AUDPC, col=AUDPC))+
    geom_tile()+
    facet_wrap(.~genotype)+
    scale_fill_gradient2(low="#4B375C",mid = "#EED78C", high= "#FFB833", midpoint = 250)+
    scale_color_gradient2(low="#4B375C",mid = "#EED78C", high= "#FFB833", midpoint = 250)+
    theme_bw()+
    xlab("lag phase [min]")+
    ylab("LDT")+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11,  color="black",
                                     #face="italic"
          ),
          axis.title = element_text(size=11, color="Black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )
  
)
