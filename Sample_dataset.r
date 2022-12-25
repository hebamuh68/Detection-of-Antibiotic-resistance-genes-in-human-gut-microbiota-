# Read tap separated data
# human intestinal metagenome drug resistance gene study 53 genomics
# samples collected from healthy volunteers, didn't take any antibiotic agent in the past 3 months
# Signal intensities of perfect match (PM) and mismatch (MM) probes

# Age 3f, 5m
pre_school_female = read.table("data/pre school/female/GSM1306735_426998A08_20101019_0780_532.pair", sep="\t", header=TRUE)
pre_school_male = read.table("data/pre school/male/GSM1306734_426998A07_20101019_0780_532.pair", sep="\t", header=TRUE)

# Age 30f, 33m
adult_female = read.table("data/adult group/female/GSM1306672_388797A03_100121_532_532.pair", sep="\t", header=TRUE)
adult_male = read.table("data/adult group/male/GSM1306673_388797A04_100121_532_532.pair", sep="\t", header=TRUE)

# Age 16f, 18m
high_school_female = read.table("data/high school group/female/GSM1306757_453839A12_01030_0790_532.pair", sep="\t", header=TRUE)
high_school_male = read.table("data/high school group/male/GSM1306753_453839A09_01030_0790_532.pair", sep="\t", header=TRUE)

# Age 10fm
school_female = read.table("data/school group/female/GSM1306774_453840A04_110318_532.pair", sep="\t", header=TRUE)
school_male = read.table("data/school group/male/GSM1306772_453840A02_110318_532.pair", sep="\t", header=TRUE)


#===============================================log2 transform=======================================
# log2 transform
log2_transform = function(data){
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { data[which(data <= 0)] <- NaN
  new_data <- log2(data) }
  return(new_data)
}

pre_school = c(pre_school_female$PM, pre_school_male$PM)
Npre_school = log2_transform(pre_school)
pre_school_mean = mean(Npre_school)


school = c(school_female$PM, school_male$PM)
Nschool = log2_transform(school)
school_T = t.test(Nschool, mu = pre_school_mean) 
school_T

high_school = c(high_school_female$PM, high_school_male$PM)
Nhigh_school = log2_transform(high_school)
high_school_mean = mean(high_school)
high_school_T = t.test(Nhigh_school, mu = pre_school_mean) 
high_school_T

adult = c(adult_female$PM, adult_male$PM)
Nadult = log2_transform(adult)
adult_mean = mean(adult)
adult_T = t.test(Nadult, mu = pre_school_mean) 
adult_T

data_matrix = matrix(ncol=4)
data_matrix = cbind(Npre_school, Nschool, Nhigh_school, Nadult)

#====================================box-and-whisker plot============================================
png(filename="sample_plots/Sample_box-and-whisker.png",width=850)
par(mar=c(7,4,2,1))
title <- paste ("GSE54070", "/", annotation(gset), sep ="")
boxplot(data_matrix, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
dev.off()

#====================================expression value distribution plot=============================
png(filename="sample_plots/Sample_expression_value.png",width=850)
par(mar=c(4,4,2,1))
plotDensities(data_matrix,group=NULL, main=title, col=c("red","blue", "green", "black"))
dev.off()

#=================================== mean-variance trend =============================
png(filename="sample_plots/mean-variance.png",width=850)
fit = lmFit(data_matrix)
plotSA(fit, main="Mean variance trend, GSE54070")
dev.off()
