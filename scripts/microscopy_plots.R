library(ggplot2)
library(reshape2)
library(dplyr)
library(rstatix)

set.seed(15091998)
setwd("PHD/Projets/T7/Figures/Ascet/")

lateral <- read.csv(file = 'lateral.csv', sep = ';')
cols <- c("CFP_T7_OFF", "CFP_T7_ON", "CFP_T7_OFF_RIF", "CFP_T7_ON_RIF", "YFP_T7_OFF", "YFP_T7_ON", "YFP_T7_OFF_RIF", "YFP_T7_ON_RIF")
colnames(lateral) <- cols

lateral_stack <- stack(lateral)

# Do t-test on all of them
group_by(lateral_stack, ind) %>%
  summarise(median(values, na.rm = TRUE))
stat.test <- lateral_stack %>%
  t_test(values ~ ind) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()


ggplot(lateral_stack, aes(x=values, color=ind)) + 
  geom_density(na.rm = TRUE)


longitidunal_1focus <- read.csv(file = 'longitudinal_1focus.csv', sep = ';')
colnames(longitidunal_1focus) <- c("CFP_T7_OFF", "CFP_T7_ON", "CFP_T7_OFF_RIF", "CFP_T7_ON_RIF", "YFP_T7_OFF", "YFP_T7_ON", "YFP_T7_OFF_RIF", "YFP_T7_ON_RIF")

longitidunal_1focus_stack <- stack(longitidunal_1focus)
longitidunal_1focus_stack$values_rev <- 1 - longitidunal_1focus_stack$values
longitidunal_1focus_stack <- transform(longitidunal_1focus_stack, values_corr=pmin(values, values_rev))
longitidunal_1focus_stack$values_corr[longitidunal_1focus_stack$values_corr < 0] <- NA
ggplot(longitidunal_1focus_stack, aes(x=values_corr, color=ind)) + 
  geom_density(na.rm = TRUE)

# Do t-test on all of them
group_by(longitidunal_1focus_stack, ind) %>%
  summarise(median(values_corr, na.rm = TRUE))
stat.test <- longitidunal_1focus_stack %>%
  t_test(values_corr ~ ind) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

longitidunal_2focus <- read.csv(file = 'longitudinal_2foci.csv', sep = ';')
colnames(longitidunal_2focus) <- c("CFP_T7_OFF", "CFP_T7_ON", "CFP_T7_OFF_RIF", "CFP_T7_ON_RIF", "YFP_T7_OFF", "YFP_T7_ON", "YFP_T7_OFF_RIF", "YFP_T7_ON_RIF")

longitidunal_2focus_stack <- stack(longitidunal_2focus)
longitidunal_2focus_stack$values[longitidunal_2focus_stack$values < 0] <- NA
longitidunal_2focus_stack$values[longitidunal_2focus_stack$values > 1] <- NA
longitidunal_2focus_stack$values_rev <- 1 - longitidunal_2focus_stack$values
longitidunal_2focus_stack <- transform(longitidunal_2focus_stack, values_corr=pmin(values, values_rev))
longitidunal_2focus_stack$values_corr[longitidunal_2focus_stack$values_corr < 0] <- NA

# Do t-test on all of them
group_by(longitidunal_2focus_stack, ind) %>%
  summarise(median(values_corr, na.rm = TRUE))
stat.test <- longitidunal_2focus_stack %>%
  t_test(values_corr ~ ind) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

ggplot(longitidunal_2focus_stack, aes(x=values, color=ind)) + 
  geom_density(na.rm = TRUE)



for (col in cols) {
  #Generate Simple Kernel Density Estimate uing the default R function
  fit1 <- density(x=c(lateral[col]), na.rm = TRUE)    # Replace by the curve  
  
  #Bootstrap starts
  fit2 <- replicate(1000,
                    { 
                      #Sample with replacement, for the bootstrap from the
                      #original dataset and save the resample to x
                      x <- sample(c(lateral[col]), replace=TRUE)                    
                      
                      #Generate the density from the resampled dataset, and
                      #extract y coordinates to generate variablity bands
                      #for that particular x coordinate in the smooth curve
                      density(x, from=min(fit1$x), to=max(fit1$x), na.rm = TRUE)$y
                    }) 
  
  #Apply the quantile function to the y coordinates to get the
  #bounds of the polygon to be drawn on the y axis
  #if so, why the 2.5% to 97.5% range
  fit3 <- apply(fit2, 1, quantile, c(0.025,0.975) )
  
  matplot(fit1$x, cbind(fit1$y, t(fit3)), type = "l",
          lty = c(1, 2, 2), lwd = 3, col = c("black", "red", "red"))
  plot(fit1, ylim = range(fit3))
  polygon( c(fit1$x, rev(fit1$x)), c(fit3[1,], rev(fit3[2,])),
           col='black', density = -0.5, border=F)
  lines(fit1, col = "red", lwd = 3)
}

longitidunal_1focus$values_rev <- 1 - longitidunal_1focus$YFP_T7_ON_RIF
longitidunal_1focus <- transform(longitidunal_1focus, values_corr=pmin(YFP_T7_ON_RIF, values_rev))
longitidunal_1focus$values_corr[longitidunal_1focus$values_corr < 0] <- NA

#Generate Simple Kernel Density Estimate uing the default R function
fit1 <- density(x=c(longitidunal_1focus$values_corr), na.rm = TRUE)    # Replace by the curve  

#Bootstrap starts
fit2 <- replicate(1000,
                  { 
                    #Sample with replacement, for the bootstrap from the
                    #original dataset and save the resample to x
                    x <- sample(c(longitidunal_1focus$values_corr), replace=TRUE)                    
                    
                    #Generate the density from the resampled dataset, and
                    #extract y coordinates to generate variablity bands
                    #for that particular x coordinate in the smooth curve
                    density(x, from=min(fit1$x), to=max(fit1$x), na.rm = TRUE)$y
                  }) 

#Apply the quantile function to the y coordinates to get the
#bounds of the polygon to be drawn on the y axis
#if so, why the 2.5% to 97.5% range
fit3 <- apply(fit2, 1, quantile, c(0.025,0.975) )

matplot(fit1$x, cbind(fit1$y, t(fit3)), type = "l",
        lty = c(1, 2, 2), lwd = 3, col = c("black", "red", "red"))
plot(fit1, ylim = range(fit3))
polygon( c(fit1$x, rev(fit1$x)), c(fit3[1,], rev(fit3[2,])),
         col='black', density = -0.5, border=F)
lines(fit1, col = "red", lwd = 3)



ks.test(lateral$CFP_T7_OFF, lateral$YFP_T7_OFF)
ks.test(longitidunal_2focus$CFP_T7_OFF, longitidunal_2focus$YFP_T7_OFF)
ks.test(longitidunal_1focus$CFP_T7_OFF, longitidunal_1focus$YFP_T7_OFF)
ks.test(lateral$CFP_T7_ON, lateral$YFP_T7_ON)
ks.test(lateral$YFP_T7_OFF, lateral$YFP_T7_OFF_RIF)
ks.test(lateral$CFP_T7_ON, lateral$CFP_T7_OFF)
ks.test(longitidunal_2focus$YFP_T7_OFF, longitidunal_2focus$YFP_T7_OFF_RIF)
ks.test(longitidunal_2focus$YFP_T7_OFF, longitidunal_2focus$YFP_T7_ON_RIF)
ks.test(longitidunal_1focus$YFP_T7_OFF, longitidunal_1focus$YFP_T7_OFF_RIF)
ks.test(longitidunal_1focus$YFP_T7_OFF, longitidunal_1focus$YFP_T7_ON_RIF)

interfocal_distance <- read.csv('interfocal_distance.csv', sep=";")
t.test(interfocal_distance$ARA, interfocal_distance$NI) 
t.test(interfocal_distance$ARA, interfocal_distance$ARA.RIF) 
t.test(interfocal_distance$ARA, interfocal_distance$RIF) 
t.test(interfocal_distance$NI, interfocal_distance$ARA.RIF) 
t.test(interfocal_distance$NI, interfocal_distance$RIF) 
t.test(interfocal_distance$ARA.RIF, interfocal_distance$RIF)

interfocal_distance_stack = stack(interfocal_distance)
interfocal_distance_stack = mutate(interfocal_distance_stack, ind = factor(ind, levels=c("NI", "ARA", "RIF", "ARA.RIF")))

ggplot(interfocal_distance_stack, aes(y=values, x=ind, color=ind)) +
  geom_point() + 
  ylim(0.75, 1) +
  theme_minimal()
summary(interfocal_distance)
