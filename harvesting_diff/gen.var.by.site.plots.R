##### R script to create lm plots on pi, He, MAHF, gen load
setwd("~/Desktop/PD_stuffies/R_work/gen.div")
library(ggplot2)
library(cowplot)

## read in data
lat_by_gen <- read_excel("lat.by.gen.xlsx", col_types = c("text", "numeric", "numeric",  "numeric", "numeric"))

#gen.df<- lat_by_gen[rev(rownames(lat_by_gen)),]

#get site order
lat_by_gen$Site <- factor(lat_by_gen$Site , levels=c("KR", "FR", "BB", "DB", "HP", "CR", "SB", "VN", "GP", "CA", "WA", "IP", "SH", "SC", "CB", "PB", "SR", "BT", "BA"), ordered = TRUE)
#gen.df$Site=factor(gen.df$Site, levels=c("BA", "BT", "SR", "PB", "CB", "SC", "SH", "IP", "WA", "CA", "GP", "VN", "SB", "CR", "HP", "DB", "BB", "FR", "KR"), ordered = TRUE)

#name colors
#ALL.cols=c("#436EEE", "#6959CD","#009ACD","#B9D3EE","#40E0D0","#458B00","#66CD00", "#00FF00","#8B7500","#FFFF00", "#FFD700", "#CD3700", "#FF1493","#FF83FA", "#8B1C62", "#CD6090", "#9F79EE", "#68228B")  

ALL.cols=c("#436EEE", "#6959CD","#009ACD","#B9D3EE","#40E0D0","#458B00","#66CD00", "#00FF00","#8B7500","#FFFF00", "#FFD700","#EE7600", "#CD3700", "#FF1493","#FF83FA", "#8B1C62", "#CD6090", "#9F79EE", "#68228B")

#create plots per variable
pi <- ggplot(lat_by_gen, aes(x=Lat, y=pi)) + geom_point(aes(color=Site), size=5)+geom_point(shape = 1,size = 5,colour = "black") + geom_smooth(color="black", fullrange=TRUE) +scale_color_manual(values=ALL.cols)+theme_minimal()+theme(axis.text=element_text(size=14))

He <- ggplot(lat_by_gen, aes(x=Lat, y=He)) + geom_point(aes(color=Site),  size=5)+geom_point(shape = 1,size = 5,colour = "black") + geom_smooth(color="black", fullrange=TRUE) +scale_color_manual(values=ALL.cols)+theme_minimal()+theme(axis.text=element_text(size=14))

D <- ggplot(lat_by_gen, aes(x=Lat, y=D)) + geom_point(aes(color=Site),  size=5)+geom_point(shape = 1,size = 5,colour = "black") + geom_smooth(color="black", fullrange=TRUE) +scale_color_manual(values=ALL.cols)+theme_minimal()+theme(axis.text=element_text(size=14))

theta <- ggplot(lat_by_gen, aes(x=Lat, y=theta)) + geom_point(aes(color=Site),  size=5)+geom_point(shape = 1,size = 5,colour = "black") + geom_smooth(color="black", fullrange=TRUE) +scale_color_manual(values=ALL.cols)+theme_minimal()+theme(axis.text=element_text(size=14))

MAHF<- ggplot(lat_by_gen, aes(x=Lat, y=MAHF)) + geom_point(aes(color=Site),  size=5)+geom_point(shape = 1,size = 5,colour = "black") + geom_smooth(method=lm, color="black", fullrange=TRUE) +scale_color_manual(values=ALL.cols)+theme_minimal()+theme(axis.text=element_text(size=14))

gen_load <- ggplot(lat_by_gen, aes(x=Lat, y=gen_load)) + geom_point(aes(color=Site),  size=5)+geom_point(shape = 1,size = 5,colour = "black") + geom_smooth(color="black", fullrange=TRUE) +scale_color_manual(values=ALL.cols)+theme_minimal()+theme(axis.text=element_text(size=14))

prv_snps <- ggplot(lat_by_gen, aes(x=Lat, y=prv_snps)) + geom_point(aes(color=Site),  size=5)+geom_point(shape = 1,size = 5,colour = "black") + geom_smooth(method=lm, color="black", fullrange=TRUE) +scale_color_manual(values=ALL.cols)+theme_minimal()+theme(axis.text=element_text(size=14))

#get legend from a plot
legend <- get_legend(pi)

#arrange plots (for MS fig)
pgrd<-plot_grid(pi+theme(legend.position="none"), gen_load+theme(legend.position="none"), prv_snps+theme(legend.position="none"), MAHF+theme(legend.position="none"), labels = c('A', 'B','C', 'D'))


#arrange plots (for SM fig)
pgrd<-plot_grid(D+theme(legend.position="none"), theta+theme(legend.position="none"), labels = c('A', 'B'))


pdf("gen.div.NEW.NEW.pdf")
plot_grid(pgrd, legend, rel_widths = c(3, .4))
dev.off()

pdf("pi.GAM.scatter.plots.pdf")
pi
dev.off()

pdf("He.GAM.scatter.plots.pdf")
He
dev.off()

pdf("D.GAM.scatter.plots.pdf")
D
dev.off()

pdf("theta.GAM.scatter.plots.pdf")
theta
dev.off()

pdf("MAHF.GAM.scatter.plots.pdf")
MAHF
dev.off()


pdf("gen_load.GAM.scatter.plots.pdf")
gen_load
dev.off()


## do linear plot for MAHF 
MAHF.lm <- ggplot(lat_by_gen, aes(x=Lat, y=MAHF)) +geom_point(aes(color=Site), size=3)+ geom_smooth(method=lm, color="black", fullrange=TRUE) + scale_color_manual(values=ALL.cols)+theme_minimal()+scale_x_reverse()+theme(axis.text=element_text(size=14))

pdf("MAHF.LM.scatter.plots.pdf")
MAHF.lm
dev.off()

## and for prv snps
prv.snp.lm <- ggplot(lat_by_gen, aes(x=Lat, y=prv_snps)) +geom_point(aes(color=Site), size=3)+ geom_smooth(method=lm, color="black", fullrange=TRUE) + scale_color_manual(values=ALL.cols)+theme_minimal()+scale_x_reverse()+theme(axis.text=element_text(size=14))
pdf("prv_snps.LM.scatter.plots.pdf")
prv.snp.lm
dev.off()

#run linear models per variable
fit.pi <- lm(pi ~ Lat, data = lat_by_gen)
summary(fit.pi)


Coefficients:
  Estimate Std. Error t value Pr(>|t|)   
(Intercept)  0.171382   0.049045   3.494  0.00278 **
  Lat         -0.000729   0.001433  -0.509  0.61755   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.02196 on 17 degrees of freedom
Multiple R-squared:  0.01499,	Adjusted R-squared:  -0.04295 
F-statistic: 0.2587 on 1 and 17 DF,  p-value: 0.6176

##
fit.theta <- lm(theta ~ Lat, data = lat_by_gen)
summary(fit.theta)


Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.157e-03  7.323e-05  15.798 1.36e-11 ***
  Lat         -2.049e-06  2.140e-06  -0.957    0.352    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 3.279e-05 on 17 degrees of freedom
Multiple R-squared:  0.05116,	Adjusted R-squared:  -0.004649 
F-statistic: 0.9167 on 1 and 17 DF,  p-value: 0.3518

##
fit.mahf <- lm(MAHF ~ Lat, data = lat_by_gen)
summary(fit.mahf)

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept) 118399.9    17750.6   6.670 3.95e-06 ***
  Lat          -2390.1      518.8  -4.607 0.000251 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 7948 on 17 degrees of freedom
Multiple R-squared:  0.5553,	Adjusted R-squared:  0.5292 
F-statistic: 21.23 on 1 and 17 DF,  p-value: 0.000251


##
fit.gl <- lm(gen_load ~ Lat, data = lat_by_gen)
summary(fit.gl)


Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)  0.836694   0.035056  23.867 1.64e-14 ***
  Lat         -0.001708   0.001025  -1.668    0.114    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0157 on 17 degrees of freedom
Multiple R-squared:  0.1406,	Adjusted R-squared:  0.09003 
F-statistic: 2.781 on 1 and 17 DF,  p-value: 0.1137

#prv snps

summary(fit.prv)

Call:
  lm(formula = prv_snps ~ Lat, data = lat_by_gen)

Residuals:
  Min     1Q Median     3Q    Max 
-44150 -12338   1651   8400  40409 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)   251422      49392   5.090 0.000109 ***
  Lat            -5268       1442  -3.654 0.002140 ** 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 22070 on 16 degrees of freedom
Multiple R-squared:  0.4549,	Adjusted R-squared:  0.4208 
F-statistic: 13.35 on 1 and 16 DF,  p-value: 0.00214



