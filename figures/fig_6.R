# Influential factors of symptom-severity 
# using a mixed model to estimate factor-size and model AUDPC
library(tidyverse)
library(lsmeans)
library(nlme)
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

coef(mod_contr)

#write.csv(anova_coefficients, "C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/2024_01_30_coefficients_audpc_anova.csv")

(output_lme <- mod_contr2$coefficients$fixed )
#write.csv(output_lme, "C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/2024_01_30_coefficients_audpc_estimates_fixed.csv")


##
# get estimates of each species 

#AUDPC_la1282 = 1181.8118 - 0.1428 * lag - 151.3654 * LDT + 0.0169 * LDT * lag
#AUDPC_la1809 = 1545.3930 - 0.4204 * lag - 184.8083 * LDT + 0.0573 * LDT * lag
#AUDPC_la1941 = 1555.5933 - 0.3630 * lag - 215.5537 * LDT + 0.0508 * LDT * lag

#


# Define ranges for LDT
# hours to LDT(mins) log(x * 60)  
# (8640-1440)*0.01
ldt_range <- seq(5, 7.31, 0.0231)

# define lag-range 
lag_range <- seq(1440, 4320, 28.8)

genotypes <- c("LA1282", "LA1809", "LA1941")

# Generate empty dataframes to store modeled data 

audpc_values <- data.frame()



# Iterate over LDT and lag values per genotype

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



# Plot models as heatmap
coordinates_adj_mean <- data.frame(genotype=c("LA1282", "LA1809", "LA1941"),
                               mean_lag=c((2.37*(24*60)), (1.59*(24*60)),(2.08*(24*60))),
                               mean_ldt= c(6.03, 6.04, 6.38))

# coordinates of model-adjusted values
coordinates_1809 <- data.frame(mean_lag=1.59*(24*60), 
                               mean_ldt=6.04)

coordinates_1941 <- data.frame(mean_lag=2.08*(24*60), mean_ldt=6.38)





(p <- audpc_values %>% 
    mutate(lag=lag/(24*60)) %>% 
    ggplot(aes(x=lag, y=LDT))+
  geom_tile(aes(fill=AUDPC, color=AUDPC))+
    geom_point(data=coordinates_adj_mean, aes(x=mean_lag/(24*60), y=mean_ldt), col="white", 
               shape=4, size=2, stroke=2)+
  scale_fill_gradient2(low="#4B375C",mid = "#EED78C", high= "#FFB833", midpoint = 250)+
    scale_color_gradient2(low="#4B375C",mid = "#EED78C", high= "#FFB833", midpoint = 250)+
    theme_bw()+
    xlab("lag phase [d]")+
    ylab("LDT [h]")+
    scale_y_continuous(breaks=c(5.192957,5.886104, 6.291569,6.984716, 7.2732),labels = c(3,6, 9, 18, 24))+
    facet_wrap(.~genotype)+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11,  color="black",
                                     #face="italic"
          ),
          axis.title = element_text(size=11, color="Black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0),
                                      face="bold"),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'), 
          
          strip.background = element_rect(fill="white"), 
          strip.text = element_text(face="bold", size=13)
    )
)



ggsave("C:/Users/suaph281/Desktop/test_heatmap.png", 
       dpi=800, width=9, height=3, bg="white")

ggsave("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/figures/fig_6.png", 
       dpi=600, width=10, height=3, bg="white")




