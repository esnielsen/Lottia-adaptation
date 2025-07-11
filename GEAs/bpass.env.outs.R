# Get env-correlated baypass outliers

library(dplyr)

setwd("/ocean/projects/deb200006p/enielsen/LGwork/Outliers/03_baypass")

dB.df=read.table("ALL_2_r2.NEW.env.output_summary_betai_reg.out",h=T)

names(dB.df)[7] <- "dB"

dB.df.outs <- dB.df %>% filter(dB>20)

bpass.env.out.list = subset(dB.df.outs , select = -c(3:6,8:10) )

write.table(bpass.env.out.list, "bpass.NEW.env.out.list.txt", row.names=FALSE, sep="\t", quote = FALSE)
