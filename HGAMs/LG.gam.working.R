
library(mgcv)
library(MASS)
library(stringr)
library(gamm4)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(viridis)
library(cowplot)
library(kableExtra)
library(docxtools)
library(knitr)
library(tibble)
library(dplyr)
library(gratia)
library(latex2exp)


#load data and remove zero counts
size_count_vis_plot <- ggplot(dplyr::filter(LG_count_FINAL, count > 0),
                        aes(x=year, y=latitude, size=count))+
  facet_wrap(~ size_bin) +
  geom_point() +
  scale_size(name = "Count", range = c(0.2, 3)) +
  labs(x = "Year", y = "Latitude") +
  theme(legend.position = "bottom")

m1 <- gam(count ~ s(sst) + s(year) + s(latitude), data=LG_count_FINAL, method="REML", family="poisson")

m2 <- gam(count ~ size_bin + s(sst) + s(latitude,longitude), data=LG_count_FINAL, method="REML", family="poisson")

m3 <- gam(count ~ size_bin*sst + s(latitude,longitude), data=LG_count_FINAL, method="REML", family="poisson")

m4 <- gam(count ~ size_bin*sst + s(latitude,longitude, by=year_f), data=LG_count_FINAL, method="REML", family="poisson")


count_modG <- gam(count ~ te(year, latitude, k=c(20, 20)),
                 data=LG_count_FINAL, method="REML", family="poisson")

library(gratia)
draw(count_modG)

count_pred <- transform(LG_count_FINAL, modG = predict(count_modG, type="response"))

ggplot(count_pred, aes(x=modG, y=count)) +
  facet_wrap(~size_bin) +
  geom_point(alpha=0.1) +
  geom_abline() +
  labs(x="Predicted count", y="Observed count")

LG_count_FINAL$size_bin <- as.factor(LG_count_FINAL$size_bin)

count_modGI <- gam(count ~ size_bin +
                    te(year, latitude, k=c(20, 20), m=2) +
                    te(year, latitude, by=size_bin,
                       k=c(20, 20), m=1),
                  data=LG_count_FINAL, method="REML", family="poisson")

count_modI <- gam(count ~ size_bin + te(year, latitude, by=size_bin,
                    k=c(20, 20), m=2),
                 data=LG_count_FINAL, method="REML", family="poisson")

#plotting output from modGI

LG_count_FINAL_GI <- transform(LG_count_FINAL, modGI = predict(count_modGI, type="response"))

count_modGI_plot <- ggplot(data=LG_count_FINAL_GI, 
                           aes(x=year, y=latitude, fill=modGI, color=modGI)) + geom_point(size = 1.7)+
  geom_tile(size=0.25) +
  facet_wrap(~ size_bin, ncol=3) +
  scale_fill_viridis("count") +
  scale_color_viridis("count") +
  #scale_x_continuous(expand=c(0, 0), breaks=c(1, 26, 52)) +
  #scale_y_continuous(expand=c(0, 0), breaks=c(0, 30, 60)) +
  labs(x = "Year", y = "Latitude") +
  theme(legend.position="right")

ggplot(data=LG_count_FINAL_GI, aes(x=modGI, y=count)) +
  facet_wrap(~ size_bin, ncol=3) +
  geom_point(alpha=0.1) +
  geom_abline() +
  labs(x="Predicted count (model *GI*)", y= "Observed count")


# trying to include sst as well now
count_sst_modGI <- gam(count ~ size_bin +
                     te(year, sst, by=size_bin,
                        k=c(10, 20), m=1)+
                     te(year, latitude, by=size_bin,
                        k=c(10, 20), m=2),
                     data=LG_train,
                     family=Gamma(link="log"), method="REML",
                     drop.unused.levels=FALSE)

LG_count_FINAL_sstI <- transform(LG_count_FINAL, mod_sstI = predict(count_sst_modGI, type="response"))

count_modGI_plot <- ggplot(data=LG_count_FINAL_sstI, 
                           aes(x=year, y=latitude, fill=mod_sstI, color=mod_sstI)) + geom_point(size = 1.7)+
  geom_tile(size=0.25) +
  facet_wrap(~ size_bin, ncol=3) +
  scale_fill_viridis("count") +
  scale_color_viridis("count") +
  #scale_x_continuous(expand=c(0, 0), breaks=c(1, 26, 52)) +
  #scale_y_continuous(expand=c(0, 0), breaks=c(0, 30, 60)) +
  labs(x = "Year", y = "Latitude") +
  theme(legend.position="right")


count_sst_yr_lat_modI <- gam(count ~ size_bin +
                         te(year, sst, latitude, by=size_bin,
                            k=c(10, 50, 50), m=2),
                       data=LG_train,
                       family=Gamma(link="log"), method="REML",
                       drop.unused.levels=FALSE)


AIC_table <- AIC(count_modG, count_modGI, count_modI)%>%
  rownames_to_column(var= "Model")%>%
  #mutate(data_source = rep(c("CO2","bird_data"), each =5))%>%
  #group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  #dplyr::select(-data_source)%>%
  #mutate_at(.vars = vars(df,AIC, deltaAIC), 
    #        .funs = funs(round,.args = list(digits=0)))

  
#set training and test data, choosing either odd or even years
  
  
LG_count_FINAL$year_f <- as.factor(LG_count_FINAL$year)

LG_train <- subset(LG_count_FINAL, year%%2==0)
LG_test  <- subset(LG_count_FINAL, year%%2==1) 


#This function calculates the deviance of out-of-sample data,
#conditional on their mean predicted value from the model
get_deviance <- function(model, y_pred, y_obs, weights = NULL){
  stopifnot(length(y_obs)==length(y_pred))
  #We don't use the weights term in this paper, but it can be useful if
  #how well the model matters more for some sample points than others
  if(is.null(weights)) weights = rep(1, times= length(y_obs))
  #this uses the deviance residual function from the model family to
  #calculate deviances for individual points
  dev_residuals = model$family$dev.resids(y_obs, y_pred, weights)
  return(sum(dev_residuals))
}

## may potentially want to use bs="cr" for year
count_sst_yearf_G_modI <- gam(count ~ s(year, by=size_bin, k=10) + 
                           s(size_bin, bs="re") + s(size_bin, year_f, bs="re")+
                           te(sst, latitude, by=size_bin, k=c(20, 20), m=1),
                           data=LG_train,
                           family=Gamma(link="log"), method="REML",
                           drop.unused.levels=FALSE)




#### On use of te() or s()

In s() option isotropy is assumed, i.e. the same amount of smoothing is used in both directions (time and month). 

This could be reasonable for spatial fitting, or for interactions where both variables are in the same unit, but not in our case. 

For interactions between variables that should not be smoothed with the same amount, we can use tensor products (te)


##### On adding "ts" to each variable to see if they improve model fit:
Instead there are some specialised basis functions that allow shrinkage: ’ts’ and ’cs’

model9 <- gamm(logTot.P ~ s(datenr, bs='ts') + s(Month, bs='ts') + s(logRunoff, bs=’ts’) + s(Abs._F, bs='ts'), data=river)

If we add a variable that is not improving the model fit the smooth is put to a straight line && in summary will have 0 F (estimated degrees of freedom).



vis.gam(count_sst_modGI, view=c("year", "sst"), 
        cond=list(size_bin="small"),
        plot.type='contour', color='topo', main='small')

### Can try to do the plot below, with year as factor, doing before blob, and after blob
vis.gam(m1, view=c("sst", "latitude"), 
        cond=list(size_bin="small", year_f='before_B'),
        plot.type='contour', color='topo', main='small-before')

plot_smooth(count_sst_modGI, view="year", cond=list(size_bin="small"), rm.ranef=TRUE, rug=FALSE, col="red")
plot_smooth(count_sst_modGI, view="year", cond=list(size_bin="med"), rm.ranef=TRUE, rug=FALSE, col="green", add=TRUE)
plot_smooth(count_sst_modGI, view="year", cond=list(size_bin="large"), rm.ranef=TRUE, rug=FALSE, col="blue", add=TRUE)


#### OCT 2021 working with currents #####


# combining uo & vo because they look relatively correlated per year

LG_count_plot_mpa_F$size_bin <- as.factor(LG_count_plot_mpa_F$size_bin)

m2 <- gam(count ~ s(year, by=size_bin, bs="ts") + s(sst, bs="ts", by=size_bin) + s(uo, vo, bs="ts", by=size_bin) + offset(log(plots_sampled)), data=LG_count_plot_mpa_F, method="REML", family=Gamma(link="log"))

m1 <- gam(count ~ s(sst, bs="ts") + s(uo, vo, bs="ts") + s(year, bs="ts") + s(latitude, longitude, bs="ts") + offset(log(plots_sampled)), data=LG_count_plot_mpa_F, method="REML", family=nb(link="log"))


p <-ggplot(LG_count_plot_mpa_F) +
  aes(x = year, y = count, color = latitude) +
  geom_point() +
  stat_smooth() +
  theme_bw()
p + facet_wrap(~size_bin)

yr.sst_modGI <- bam(count ~ size_bin +
                     te(year, sst, k=c(20, 20), m=2) +
                     te(year, sst, by=size_bin,
                        k=c(20, 20), m=1)
                   + offset(log(plots_sampled)),
                   data=LG_count_plot_mpa_F, family=nb(link="log"), method="REML")

yr.sst_modI <- bam(count ~ size_bin + te(year, sst, by=size_bin,
                                        k=c(20, 20), m=2)
                  + offset(log(plots_sampled)),
                  data=LG_count_plot_mpa_F, family=Gamma(link="log"), method="REML")

#plotting output from modGI

LG_count_FINAL_GI <- transform(LG_count_FINAL, modGI = predict(count_modGI, type="response"))

count_modGI_plot <- ggplot(data=LG_count_FINAL_GI, 
                           aes(x=year, y=latitude, fill=modGI, color=modGI)) + geom_point(size = 1.7)+
  geom_tile(size=0.25) +
  facet_wrap(~ size_bin, ncol=3) +
  scale_fill_viridis("count") +
  scale_color_viridis("count") +
  #scale_x_continuous(expand=c(0, 0), breaks=c(1, 26, 52)) +
  #scale_y_continuous(expand=c(0, 0), breaks=c(0, 30, 60)) +
  labs(x = "Year", y = "Latitude") +
  theme(legend.position="right")

ggplot(data=LG_count_FINAL_GI, aes(x=modGI, y=count)) +
  facet_wrap(~ size_bin, ncol=3) +
  geom_point(alpha=0.1) +
  geom_abline() +
  labs(x="Predicted count (model *GI*)", y= "Observed count")


# trying to include sst as well now
yr_curr_sst_modI <- bam(count ~ size_bin +
                         te(year, sst, by=size_bin,
                             m=1)+
                         te(year, uo, vo, by=size_bin, d=c(1,2),
                            bs=c("cr", "tp"))
                       + offset(log(plots_sampled)),
                       data=LG_train,
                       family=Gamma(link="log"), method="REML",
                       drop.unused.levels=FALSE)


# trying to include interaction between size & protection
LG_count_plot_mpa_F$mpa_designation <- as.factor(LG_count_plot_mpa_F$mpa_designation)

fac.model <- bam(
  count ~
    interaction(size_bin, mpa_designation) +
    s(year, bs="tp", k=20) +
    s(year, by = size_bin, bs="tp", k=20) +
    s(year, by = mpa_designation, bs="tp", k=20) +
    s(year, by = interaction(size_bin, mpa_designation), bs="tp", k=20) + offset(log(plots_sampled)),
  data=LG_count_plot_mpa_F,
  family=nb(link="log"), method="REML"
)


fac.model.2 <- bam(
  count ~
size_bin + mpa_designation + size_bin:mpa_designation + 
  s(year, by = size_bin, k=25, bs="tp") +
  s(year, by = mpa_designation, k=25, bs="tp") +
  s(year, by = interaction(size_bin, mpa_designation), k=25, bs="tp") + offset(log(plots_sampled)),
 data=LG_count_plot_mpa_F,
  family=nb(link="log"), method="REML"
)

## Fac model including 2020
sc_bin_fac_models <- read_excel("sc.bin.fac.mod.filter.xlsx")
sc_bin_fac_models$protection <- as.factor(sc_bin_fac_models$protection)
sc_bin_fac_models$size_bin <- as.factor(sc_bin_fac_models$size_bin)

fac.model.cr.tw.k20 <- bam(count ~ interaction(size_bin, protection) + s(year, bs = "cr",k = 20) + s(year, by = size_bin, bs = "cr", k = 20) + s(year, by = protection, bs = "cr", k = 20) + s(year, by = interaction(size_bin, protection), bs = "cr", k = 20) + offset(log(plots_sampled)), data=sc_bin_fac_models, family=tw, method="REML")

# To plot fac.model!
c <- getViz(fac.model.cr.tw.k20)
print(plot(c, allTerms = TRUE), pages = 13)

# To test models (basically gives same information as gam.check():
anova(m1)

# To get over fact that time series data breaks assumptions of independence (pg 172)
correlation = corAR1 (form = ∼Time | ID) -- ### using gamm instead of gam!

# As X is a count, a Poisson, negative binomial or geometric distribution is appropriate. (p.259 pdf/ 243 actual)


pairs(LG_count_plot_mpa_F[,3:4,9:11], upper.panel = panel.smooth,
      lower.panel = panel.cor)

## Did the following to see if I need to include autocorrelation between years:
> fm1 <- gls(count ~ sst +year, LG_count_plot_mpa_F, correlation = corAR1(form = ~ 1 | year))
> summary(fm1)

Generalized least squares fit by REML
Model: count ~ sst + year 
Data: LG_count_plot_mpa_F 
AIC     BIC    logLik
35483.77 35513.7 -17736.89

Correlation Structure: AR(1)
Formula: ~1 | year 
Parameter estimate(s):
  Phi 
0.2699168 #### this is the important #, values separated by a year have correlation of 0.27 

Coefficients:
  Value Std.Error   t-value p-value
(Intercept) 1678.2565  757.0404  2.216865  0.0267
sst            7.4848    1.2034  6.219913  0.0000
year          -0.8424    0.3778 -2.229729  0.0258

Correlation: 
  (Intr) sst   
sst   0.083       
year -1.000 -0.108

Standardized residuals:
  Min         Q1        Med         Q3        Max 
-1.2443039 -0.6871613 -0.2931209  0.3707752  6.9978888 

Residual standard error: 104.0436 
Degrees of freedom: 2944 total; 2941 residual



### TRYING AGAIN
count_ti_re_model.nb = gam(count ~ s(uo, vo) + s(year) + s(sst) +
                        ti(sst, year) +
                          ti(uo, vo, year, d=c(2,1)) + offset(log(plots_sampled)),
                        data=LG_count_plot_mpa_F, family=nb(link="log"), method="REML")

count_ti_re_model.p = gam(count ~ s(uo, vo) + s(year) + s(sst) +
                             ti(sst, year) +
                             ti(uo, vo, year, d=c(2,1)) + offset(log(plots_sampled)),
                           data=LG_train, family=poisson, method="REML")

### Nb is better...
Resid. Df Resid. Dev   Df Deviance
1    1297.6      15052              
2    1219.9      71705 77.7   -56653
> AIC(count_ti_re_model.nb, count_ti_re_model.p)
df      AIC
count_ti_re_model.nb  69.69167 15191.34
count_ti_re_model.p  159.01642 80246.83

count_ti_re_model.xy = gam(count ~ s(uo, vo) + s(year) + s(longitude,latitude) +
                             ti(longitude,latitude, year, d=c(2,1)) +
                             ti(uo, vo, year, d=c(2,1)) + offset(log(plots_sampled)),
                           data=LG_count_plot_mpa_F, family=nb(link="log"), method="REML")

## xy does better than SST
> AIC(count_ti_re_model.nb, count_ti_re_model.xy)
df      AIC
count_ti_re_model.nb  81.38388 32342.27
count_ti_re_model.xy 122.36773 31779.45


count_ti_re_model.xy = gam(count ~ size_bin + s(uo, vo) + s(year) + s(longitude,latitude) +
                             ti(longitude,latitude, year, d=c(2,1)) +
                             ti(uo, vo, year, d=c(2,1)) + offset(log(plots_sampled)),
                           data=LG_count_plot_mpa_F, family=nb(link="log"), method="REML")


#te worse than ti
count_ti_re_model.xy.te = gam(count ~ size_bin + s(uo, vo) + s(year) + s(longitude,latitude) +
                             te(longitude,latitude, year, d=c(2,1)) +
                             te(uo, vo, year, d=c(2,1)) + offset(log(plots_sampled)),
                           data=LG_count_plot_mpa_F, family=nb(link="log"), method="REML")

count_ti_re_model.xy.GI = gam(count ~ size_bin +
                                te(latitude, longitude, year, d=c(2,1)) +
                                te(uo, vo, year, d=c(2,1)) +
                                te(latitude, longitude, year, d=c(2,1), by=size_bin) +
                                te(uo, vo, year, d=c(2,1), by=size_bin) +
                                offset(log(plots_sampled)),
                              data=LG_count_plot_mpa_F, family=nb(link="log"), method="REML")

count_ti_re_model.xy.I = gam(count ~ size_bin +
                                te(latitude, longitude, year, d=c(2,1), by=size_bin) +
                                te(vo, uo, year, d=c(2,1), by=size_bin) +
                                offset(log(plots_sampled)),
                              data=LG_count_plot_mpa_F, family=nb(link="log"), method="REML")
