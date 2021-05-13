# Author: Gordon Burtch
# Date: May 13, 2021
# Subject: simulating exaggeration of estimates with underpowered tests.

library(Superpower)
library(patchwork)
library(ggthemes)

# The other thing I play with here is seeing how the AMIP from 
# Broderick, T., Giordano, R., & Meager, R. (2020). An Automatic Finite-Sample Robustness Metric: Can Dropping a Little Data Change Conclusions?. arXiv preprint arXiv:2011.14999.
# Associates with a lack of statistical power. Seeing if a result is more 'fragile' when it was obtained under low power. 
library(zaminfluence)
# Define function to recover AMIPs
getAMIP <- function(mod){
  reg_infl <- ComputeModelInfluence(mod)
  grad_df <- GetTargetRegressorGrads(reg_infl, "Treat")
  influence_dfs <- SortAndAccumulate(grad_df)
  GetRegressionTargetChange(influence_dfs, "prop_removed")
}

# Setup
set.seed(1001)
N = 100000

# Generate some 'population' data.
Treat <- rbinom(size=1,n=N,p=0.5)
Y <- 1*Treat + rnorm(n=N,sd=5)
summary(lm(Y~Treat))

# Here's a possible design to test for this effect. 
design_result <- ANOVA_design(design = "2b",
                              n = 2500, 
                              mu = c(0, 1), 
                              sd = sd(Y), 
                              labelnames = c("Treat", "No", "Yes"))

# ANOVA_power() tells us our effective power for a t-test would be ... 
ANOVA_power(design_result)

# More generally, we can assess power at varying sample sizes per condition, and we see that we would need X subjects
# per condition to achieve power of 0.80 with alpha = 0.05. 
plot_power(design_result, max_n = 500, desired_power=80, alpha_level=0.05)

# WHAT MIGHT WE CONCLUDE FROM "SIGNIFICANT RESULTS" IN STUDIES LACKING SUFFICIENT POWER? 
# Let's take draws that provide power of about 0.25, 0.50, 0.90 and ~1.00.
pop <- data.frame(Y,Treat)
small <- 80*2
medium <- 200*2
sufficient <- 500*2
large <- 5000

sim_results <- data.frame(pwr = rep("",40000),est = rep(0,40000),pval = rep(0,40000), sens_sign = rep(0,40000), sens_sign_sig = rep(0,40000), sens_sig = rep(0,40000))

# 10000 simulations, under the four different scenarios. 
for (i in seq(1,40000,by=4)){
  
  # Estimate models.
  mod_sm <- lm(data=sample_n(pop,size=small),Y~Treat,x=TRUE,y=TRUE)
  mod_med <- lm(data=sample_n(pop,size=medium),Y~Treat,x=TRUE,y=TRUE)
  mod_suff <- lm(data=sample_n(pop,size=sufficient),Y~Treat,x=TRUE,y=TRUE)
  mod_lrg <- lm(data=sample_n(pop,size=large),Y~Treat,x=TRUE,y=TRUE)
  
  sum_small <- summary(mod_sm)
  sum_med <- summary(mod_med)
  sum_suff <- summary(mod_suff)
  sum_lrg <- summary(mod_lrg)
  
  sim_results[i,] <- list("Terrible (pwr = 0.25)",sum_small$coefficients[2,1],sum_small$coefficients[2,4],getAMIP(mod_sm)$prop_removed[1],getAMIP(mod_sm)$prop_removed[2],getAMIP(mod_sm)$prop_removed[3])
  sim_results[i+1,] <- list("Not Great (pwr = 0.50)",sum_med$coefficients[2,1],sum_med$coefficients[2,4],getAMIP(mod_med)$prop_removed[1],getAMIP(mod_sm)$prop_removed[2],getAMIP(mod_sm)$prop_removed[3])
  sim_results[i+2,] <- list("Reasonable (pwr = 0.88)",sum_suff$coefficients[2,1],sum_suff$coefficients[2,4],getAMIP(mod_suff)$prop_removed[1],getAMIP(mod_sm)$prop_removed[2],getAMIP(mod_sm)$prop_removed[3])
  sim_results[i+3,] <- list("Excellent (pwr ~ 1.00)",sum_lrg$coefficients[2,1],sum_lrg$coefficients[2,4],getAMIP(mod_lrg)$prop_removed[1],getAMIP(mod_sm)$prop_removed[2],getAMIP(mod_sm)$prop_removed[3])
}

sim_results$sig <- sim_results$pval < 0.05

sim_results <- sim_results %>% filter(pwr!="")

# As we expect, in general, yes, with a reduced sample, we have less precision and fewer significant estimates. 
table(sim_results$pwr,sim_results$sig,dnn=c("Power","Significant Result"))

# However, conditional on a result being significant, it is also more likely to be exaggerated when you lack power!
# This comes from the fact that we have less precision; only the most exaggerated/chance results come out statistically significant with such a small sample. 
sim_results$pwr <- factor(sim_results$pwr, levels = c("Excellent (pwr ~ 1.00)","Reasonable (pwr = 0.88)","Not Great (pwr = 0.50)","Terrible (pwr = 0.25)"))

# How our estimations look depending on whether they came back significant vs. not.
s_v_ns_ests <- ggplot(data=sim_results,aes(x=est,fill=factor(sig))) + 
  geom_histogram(alpha=0.4,color="black",bins=50) + 
  geom_vline(xintercept=mean(sim_results$est),color="orange",linetype="dashed",size=1) +
  geom_vline(xintercept=0,color="gray")+
  facet_wrap(.~pwr,nrow=1,ncol=4,dir="h") +
  xlab("Treatment Effect Estimate") + 
  ylab("Density") +
  scale_fill_manual(name="p-value",values=c("red","green"),labels=c("p >= 0.05","p < 0.05"))+
  theme_bw() +
  theme(text = element_text(family="Economica"), legend.text = element_text(family="Economica")) +
  NULL

# Plot them side by side.
(s_v_ns_ests) + plot_annotation(
  title = 'Exaggerated effect estimates in an underpowered study.',
  subtitle = 'As we reduce the size of our sample, we see stat. sig. estimates are increasingly like to be exaggerated.',
  caption = 'Notes: 10,000 sims; True treatment coefficient is 1; Cohen`s d ~ 0.2`; DGP: OLS of Y ~ Treat'
) 

# Low Power = More Fragile Result?
# Let's look at the relationship between power and the AMIP (sig result can be 'broken' by removing just X proportion of obs from the data...)
s_v_ns_sens <- ggplot() + 
  geom_density(data=sim_results %>% filter(sig==TRUE),aes(x=sens_sig,fill=pwr),alpha=0.25) + 
  xlab("AMIP: Proportion of Data We Must Drop to Make Significant Estimate Insignificant") + 
  #facet_wrap(.~pwr,nrow=1,ncol=4,dir="h")+
  theme_bw() +
  theme(text = element_text(family="Economica"), legend.text = element_text(family="Economica")) +
  NULL

s_v_ns_sens

