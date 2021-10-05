# libraries ------
library(aod)
library(Rcpp)
library(MASS)
library(ResourceSelection)
library(ggiraphExtra)
library(ggplot2)
library(ggpubr)
library(rms)
library(riskRegression)
library(plotly)
library(generalhoslem)
library(pROC)
library(gridExtra)
library(reshape2)
library(epiR)
library(relaimpo)
library(dominanceanalysis)

##################################################
################# Load input #####################
##################################################

#Load Library
setwd("C:/Users/epere/Dropbox/Trabajo/Andreas-Dravet-Prediction-Tool/Data-Final/For_submission_Neurology/")
datos <- as.data.frame(read.table("SCN1A_Training_Cohort.txt", header = TRUE, sep = "\t"))
s.dravet <- subset.data.frame(datos, subset = datos$Phenotype=="Dravet")
s.gesfp <- subset.data.frame(datos, subset = datos$Phenotype=="GEFS+")
datos.c <- subset.data.frame(datos, subset = datos$CADD!="NA")
datos.r <- subset.data.frame(datos, subset = datos$REVEL!="NA")
datos.m <- subset.data.frame(datos.r, subset = datos.r$CADD!="NA")
# Sample composition
samples.a <- as.data.frame(read.table("v2_Model_available_samples.txt", header = TRUE, sep = "\t"))
# Validation cohorts
datos.va <- as.data.frame(read.table("SCN1A_Validation_Cohort_1.txt", header = TRUE, sep = "\t"))
datos.vb <- as.data.frame(read.table("SCN1A_Validation_Cohort_2.txt", header = TRUE, sep = "\t"))


##################################################
##### FIGURE 2. Onset and SCN1A Genetic Score ####
################################################## 

wilcoxon.onset <- compare_means(Onset ~ Phenotype, data = datos)
wilcoxon.scn1ascore <- compare_means(SCN1A.score ~ Phenotype, data = datos)

p8 <- ggplot(data = datos, aes(x=Onset,fill=Phenotype)) +
  geom_density(alpha=0.75) +
  scale_fill_manual(labels=c("Dravet Syndrome", "GEFS+"),
                    breaks=c("Dravet","GEFS+"),
                    values=c("mediumpurple4","grey48")) +
  theme_bw() +
  theme(legend.box.background = element_rect(),
        legend.position = c(0.8, 0.4),
        text = element_text(size=8))  +
  labs(title="Age of seizure onset in Training cohort",
       x = "Age of seizure onset",
       y = "Fraction of Patients ")

p9 <- ggplot(data = datos, aes(x=SCN1A.score,fill=Phenotype)) +
  geom_density(alpha=0.75,show.legend = F) +
  scale_fill_manual(labels=c("Dravet Syndrome", "GEFS+"),
                    breaks=c("Dravet","GEFS+"),
                    values=c("mediumpurple4","grey48")) +
  theme_bw() +
  theme(text = element_text(size=8)) +
  labs(title="SCN1A genetic score in Training cohort",
       x = "SCN1A genetic score",
       y = "Fraction of Patients ")

lay <- rbind(c(1),
             c(2))
Figure2 <- grid.arrange(grobs = list(p8,p9), layout_matrix = lay)
ggsave(filename = "Figure-2.pdf", plot = Figure2, width = 3.3, height = 4.6)

##################################################
################# Main modeling ##################
##################################################

# SCN1A score & Onset model (Final) ----------------------------
# Extracting Variables
r<-as.numeric(datos$RESPONSE)
o<-as.numeric(datos$Onset)
p<-as.numeric(datos$SCN1A.score)
# Required, not clear what it does 
dd = datadist(r, o, p)
options(datadist='dd')
# Modeling
m<-lrm(r~rcs(o,3) + p, x=TRUE , y=TRUE)
s = Score(list('model' = m), 
          formula=RESPONSE~1,
          data=datos,
          plots=c("roc", "calibration"),
          B=1000, trainseeds = 1:1000, summary="IPA")
# Extracting values
modelo = s
auc = modelo$AUC$score$AUC*100
auc.ci1 = modelo$AUC$score$lower*100
auc.ci2 = modelo$AUC$score$upper*100
# IPA 
ipa = modelo$Brier$score$IPA[2]*100
brier = modelo$Brier$score$Brier[1]*100
up = modelo$Brier$contrasts$upper*100*-1
low = modelo$Brier$contrasts$lower*100*-1
ipa.ci1 = up/brier*100
ipa.ci2 = low/brier*100
stats.raw = round(cbind(auc, auc.ci1, auc.ci2, ipa, ipa.ci1, ipa.ci2),2)
stats = cbind("SCN1A genetic score + onset", stats.raw)
modelo.results = rbind(stats)


# Onset-Only Model --------------------------------------------------------
mo<-lrm(r~rcs(o,3), x=TRUE , y=TRUE)
so = Score(list('model' = mo), 
           formula=RESPONSE~1,
           data=datos,
           plots=c("roc", "calibration"),
           B=1000, trainseeds = 1:1000, summary="IPA")

# Extracting value
modelo = so
auc = modelo$AUC$score$AUC*100
auc.ci1 = modelo$AUC$score$lower*100
auc.ci2 = modelo$AUC$score$upper*100
# IPA 
ipa = modelo$Brier$score$IPA[2]*100
brier = modelo$Brier$score$Brier[1]*100
up = modelo$Brier$contrasts$upper*100*-1
low = modelo$Brier$contrasts$lower*100*-1
ipa.ci1 = up/brier*100
ipa.ci2 = low/brier*100
stats.raw = round(cbind(auc, auc.ci1, auc.ci2, ipa, ipa.ci1, ipa.ci2),2)
stats = cbind("Onset only", stats.raw)
modelo.results = rbind(modelo.results, stats)

# CADD-Onset Model --------------------------------------------------------
# Made on the subset of available data
r<-as.numeric(datos.c$RESPONSE)
o<-as.numeric(datos.c$Onset)
cadd<-as.numeric(datos.c$CADD)
# Required, not clear what it does 
dd = datadist(r, o, cadd)
options(datadist='dd')
# Modeling
mc<-lrm(r~rcs(o,3)+cadd, x=TRUE , y=TRUE)
sc = Score(list('model' = mc), 
           formula=RESPONSE~1,
           data=datos.c,
           plots=c("roc", "calibration"),
           B=1000, trainseeds = 1:1000, summary="IPA")
# Extracting value
modelo = sc
auc = modelo$AUC$score$AUC*100
auc.ci1 = modelo$AUC$score$lower*100
auc.ci2 = modelo$AUC$score$upper*100
# IPA 
ipa = modelo$Brier$score$IPA[2]*100
brier = modelo$Brier$score$Brier[1]*100
up = modelo$Brier$contrasts$upper*100*-1
low = modelo$Brier$contrasts$lower*100*-1
ipa.ci1 = up/brier*100
ipa.ci2 = low/brier*100
stats.raw = round(cbind(auc, auc.ci1, auc.ci2, ipa, ipa.ci1, ipa.ci2),2)
stats = cbind("CADD + onset", stats.raw)
modelo.results = rbind(modelo.results, stats)

# REVEL model -------------------------------------------------------------
# Made on the subset of available data
r<-as.numeric(datos.r$RESPONSE)
o<-as.numeric(datos.r$Onset)
revel<-as.numeric(datos.r$REVEL)
# Required, not clear what it does 
dd = datadist(r, o, revel)
options(datadist='dd')
# Modeling
mv<-lrm(r~rcs(o,3)+revel, x=TRUE , y=TRUE)
sv = Score(list('model' = mv), 
           formula=RESPONSE~1,
           data=datos.r,
           plots=c("roc", "calibration"),
           B=1000, trainseeds = 1:1000, summary="IPA")
# Extracting value
modelo = sv
auc = modelo$AUC$score$AUC*100
auc.ci1 = modelo$AUC$score$lower*100
auc.ci2 = modelo$AUC$score$upper*100
# IPA 
ipa = modelo$Brier$score$IPA[2]*100
brier = modelo$Brier$score$Brier[1]*100
up = modelo$Brier$contrasts$upper*100*-1
low = modelo$Brier$contrasts$lower*100*-1
ipa.ci1 = up/brier*100
ipa.ci2 = low/brier*100
stats.raw = round(cbind(auc, auc.ci1, auc.ci2, ipa, ipa.ci1, ipa.ci2),2)
stats = cbind("REVEL + onset", stats.raw)
modelo.results = rbind(modelo.results, stats)
modelo.results


##################################################
###### Figure 3 & 4(Plots-ROC-Calibration) #######
################################################## 

plotCalibration(s, pseudo = FALSE, rug = FALSE, col = "dodgerblue3")
plotCalibration(so, pseudo = FALSE, rug = FALSE, col = "darkorange1")
plotCalibration(sc, pseudo = FALSE, rug = FALSE, col = "forestgreen")
plotCalibration(sv, pseudo = FALSE, rug = FALSE, col = "purple")

# ggplot version
s.data <- as.data.frame(s$ROC$plotframe)
so.data <- as.data.frame(so$ROC$plotframe)
sc.data <- as.data.frame(sc$ROC$plotframe)
sv.data <- as.data.frame(sv$ROC$plotframe)

s.data$model <- "SCN1A genetic score + Onset"
so.data$model <- "Onset-only"
sc.data$model <- "CADD + Onset"
sv.data$model <- "REVEL + Onset"

# TPR is sensitivity(y-axis), 1-FPR is =1-Specificity (x-axis)
s.all <- as.data.frame(rbind(s.data, so.data, sc.data, sv.data))
s.all$model <- factor(s.all$model, levels = c("SCN1A genetic score + Onset", "Onset-only","CADD + Onset","REVEL + Onset"))
Figure3 <- ggplot(data = s.all, aes(x=(s.all$FPR), y=s.all$TPR,colour=s.all$model)) +
  scale_colour_manual(values=c("dodgerblue3", "darkorange1", "forestgreen", "purple", "gold")) +
  geom_line(size=0.5,show.legend = FALSE) +
  geom_segment(aes(x = 0, y = 0, xend = 0.234, yend = 0.833),size=0.5, colour = "grey44") +
  geom_segment(aes(x = 0, y = 0, xend = 0.234, yend = 0.833),size=0.5, colour = "grey44") +
  geom_segment(aes(x = 0.234, y = 0.833, xend = 1, yend = 1),size=0.5, colour = "grey44") +
  geom_abline() +
  theme_bw() +
  theme(text = element_text(size=8)) +
  labs(title="ROC curve per model",
       x = "1-Specificity",
       y = "Sensitivity",
       colour = element_blank())
ggsave(filename = "Figure-3.pdf", plot = Figure3, width = 3.3, height = 3.3)
Figure3

##################################################
####### Predictions Validation Cohorts  ########## 
##################################################

#Australia
va.onset.pg <- data.frame(o=datos.va$Onset,p=datos.va$SCN1A.score)
va.predict <- predict(m, va.onset.pg, type="fitted.ind")
va.predict.f <- sprintf("%.2f", 100*va.predict)
final.va <-data.frame("Sample_ID"=datos.va$Sample,
                      "Diagnosis"=datos.va$Phenotype,
                      "Onset"=datos.va$Onset,
                      "SCN1A genetic score"=sprintf("%.2f", datos.va$SCN1A.score),
                      "Dravet_Probability"=va.predict.f)
write.table(final.va,file="V1-Australia.prediction",sep = "\t",quote = FALSE,row.names = TRUE)
final.va.o <- final.va
final.va.o$Dravet_Probability <- as.numeric(as.character(final.va.o$Dravet_Probability))
final.va.o <- final.va.o[order(final.va.o[,5], decreasing = FALSE),]
final.va.o$Patients <- 1:203

#Belgium
vb.onset.pg <- data.frame(o=datos.vb$Onset,p=datos.vb$SCN1A.score)
vb.predict <- predict(m, vb.onset.pg, type="fitted.ind")
vb.predict.f <- sprintf("%.2f", 100*vb.predict)
final.vb <-data.frame("Sample_ID"=datos.vb$Sample,
                      "Diagnosis"=datos.vb$Phenotype,
                      "Onset"=datos.vb$Onset,
                      "SCN1A genetic score"=sprintf("%.2f", datos.vb$SCN1A.score),
                      "Dravet_Probability"=vb.predict.f)
write.table(final.vb,file="V2-Belgium.prediction",sep = "\t",quote = FALSE,row.names = TRUE)
final.vb.o <- final.vb
final.vb.o$Dravet_Probability <- as.numeric(as.character(final.vb.o$Dravet_Probability))
final.vb.o <- final.vb.o[order(final.vb.o[,5], decreasing = FALSE),]
final.vb.o$Patients <- 1:72

#Combined
combined <- rbind(datos.va,datos.vb)
combined.onset.pg <- data.frame(o=combined$Onset,p=combined$SCN1A.score)
combined.predict <- predict(m, combined.onset.pg, type="fitted.ind")
combined.response <- combined$RESPONSE

# Onset only
va.onset.onset <- data.frame(o=datos.va$Onset)
va.predict.o <- predict(mo, va.onset.onset, type="fitted.ind")

vb.onset.onset <- data.frame(o=datos.vb$Onset)
vb.predict.o <- predict(mo, vb.onset.onset, type="fitted.ind")

combined.onset <- data.frame(o=combined$Onset)
combined.predict.o <- predict(mo, combined.onset, type="fitted.ind")


##################################################
#################### FIGURE 5 ####################
##################################################
# Validation Cohort Plot BARPLOT AUSTRALIA 
p1 <-ggplot() +
  geom_col(data = final.va.o, aes(x=Patients, y=Dravet_Probability,
                                  fill=Diagnosis,
                                  position = "dodge"),width = 1.1) +
  scale_fill_manual(labels=c("Dravet Syndrome", "GEFS+"),
                    breaks=c("Dravet","GEFS+"),
                    values=c("mediumpurple4","grey48")) +
  scale_y_continuous(breaks = round(seq(0, 100, by = 10),1)) +
  theme_bw() +
  theme(legend.box.background = element_rect(),
        legend.position = c(0.895, 0.2),
        text = element_text(size=9), axis.text = element_text(size=9)) +
  labs(title="A. Model predictions for Validation Cohort 1: n=203, AUC= 94.1 [91.0-97.3], IPA= 56.6 [43.5-69.7]",
       x = "Patients",
       y = "Probability of Dravet Syndrome (%)") +
  # Coord_cartesian(expand = FALSE) +
  geom_text(aes(-5, 55, label = "Probably Dravet", vjust = 1, hjust = 0), size = 3) +
  geom_text(aes(-5, 45, label = "Probably GEFS+", vjust = 0, hjust = 0), size = 3) +
  geom_hline(aes(yintercept=(50)), colour="black", linetype="dashed", show.legend = F, size = 0.5) +
  coord_cartesian(expand = FALSE)

# Validation Cohort Plot BARPLOT BELGIUM
b1 <-ggplot() +
  geom_col(data = final.vb.o, aes(x=Patients, y=Dravet_Probability,
                                  fill=Diagnosis,
                                  position = "dodge"),width = 1.1,show.legend = FALSE) +
  scale_fill_manual(labels=c("Dravet Syndrome", "GEFS+"),
                    breaks=c("Dravet","GEFS+"),
                    values=c("mediumpurple4","grey48")) +
  scale_y_continuous(breaks = round(seq(0, 100, by = 10),1)) +
  theme_bw() +
  theme(text = element_text(size=9),axis.text = element_text(size=9)) +
  labs(title="B. Model predictions for Validation Cohort 2: n=72, AUC= 92.6 [82.7-100], IPA= 58.2 [23.6-92.9]",
       x = "Patients",
       y = "Probability of Dravet Syndrome (%)") +
  # Coord_cartesian(expand = FALSE) +
  geom_text(aes(-5, 55, label = "Probably Dravet", vjust = 1, hjust = 0), size = 3) +
  geom_text(aes(-5, 45, label = "Probably GEFS+", vjust = 0, hjust = 0), size = 3) +
  geom_hline(aes(yintercept=(50)), colour="black", linetype="dashed", show.legend = F, size = 0.5) +
  coord_cartesian(expand = FALSE)


lay3 <- rbind(c(1),c(2))
Figure5 <- grid.arrange(grobs = list(p1,b1), layout_matrix = lay3)
ggsave(filename = "Figure-5.pdf", plot = Figure5, width = 6.8, height = 6.5)


##################################################
#################### FIGURE 6 ####################
##################################################

# Validation Cohort Density Plot AUSTRALIA 
p2 <-ggplot(data = final.va.o, aes(x=Dravet_Probability,fill=Diagnosis)) +
  geom_density(alpha=0.75,show.legend = F) +
  geom_vline(mapping = aes(xintercept = 50), linetype="dashed") +
  scale_fill_manual(labels=c("Dravet Syndrome", "GEFS+"),
                    breaks=c("Dravet","GEFS+"),
                    values=c("mediumpurple4","grey48")) +
  scale_x_continuous(breaks = round(seq(0, 100, by = 10),1)) +
  theme_bw() +
  theme(text = element_text(size=8)) + 
  labs(title="Validation cohort 1",
       x = "Probability of Dravet Syndrome (%)",
       y = "Fraction of Patients ") +
  facet_wrap(. ~ final.va.o$Diagnosis, dir = "v")

# Validation Cohort Density Plot BELGIUM-
b2 <-ggplot(data = final.vb.o, aes(x=Dravet_Probability,fill=Diagnosis)) +
  geom_density(alpha=0.75,show.legend = F) +
  geom_vline(mapping = aes(xintercept = 50), linetype="dashed") +
  scale_fill_manual(labels=c("Dravet Syndrome", "GEFS+"),
                    breaks=c("Dravet","GEFS+"),
                    values=c("mediumpurple4","grey48")) +
  scale_x_continuous(breaks = round(seq(0, 100, by = 10),1)) +
  theme_bw() +
  theme(text = element_text(size=8)) +
  labs(title="Validation cohort 2",
       x = "Probability of Dravet Syndrome (%)",
       y = "Fraction of Patients ") +
  facet_wrap(. ~ final.vb.o$Diagnosis, dir = "v")

lay4 <- rbind(c(1,2))
Figure6 <- grid.arrange(grobs = list(p2,b2), layout_matrix = lay4)
ggsave(filename = "Figure-6.pdf", plot = Figure6, width = 6.8, height = 4.6)

##################################################
####################   Table 1  ##################
##################################################

tabla.1 <- c("Sensitivity","ci1","ci2","Specificity","ci1","ci2","PPV","ci1","ci2","NPV","ci1","ci2","Accuracy","ci1","ci2")
# 0.5
mroc.v1 <- roc(as.numeric(datos.va$RESPONSE), va.predict)
cont <- coords(mroc.v1, .5, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.6
cont <- coords(mroc.v1, .6, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.7
cont <- coords(mroc.v1, .7, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.8
cont <- coords(mroc.v1, .8, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.9
cont <- coords(mroc.v1, .9, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))

############ vb
# 0.5
mroc.v1 <- roc(as.numeric(datos.vb$RESPONSE), vb.predict)
cont <- coords(mroc.v1, .5, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.6
cont <- coords(mroc.v1, .6, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.7
cont <- coords(mroc.v1, .7, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.8
cont <- coords(mroc.v1, .8, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.9
cont <- coords(mroc.v1, .9, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))

################# combined
# 0.5
mroc.v1 <- roc(as.numeric(combined$RESPONSE), combined.predict)
cont <- coords(mroc.v1, .5, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.6
cont <- coords(mroc.v1, .6, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.7
cont <- coords(mroc.v1, .7, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.8
cont <- coords(mroc.v1, .8, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))
# 0.9
cont <- coords(mroc.v1, .9, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.1 <- rbind(tabla.1, sprintf("%.1f", 100*valores))

write.table(tabla.1,file = 'Table_1',quote = FALSE,sep = "\t")


##################################################
#################### Supp Table ##################
##################################################

sup.tabla.1 <- c("Sensitivity","ci1","ci2","Specificity","ci1","ci2","PPV","ci1","ci2","NPV","ci1","ci2","Accuracy","ci1","ci2")
# 0.5
mroc.v1 <- roc(as.numeric(datos.va$RESPONSE), va.predict.o)
cont <- coords(mroc.v1, .5, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.6
cont <- coords(mroc.v1, .6, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.7
cont <- coords(mroc.v1, .7, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.8
cont <- coords(mroc.v1, .8, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.9
cont <- coords(mroc.v1, .9, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))

############ vb
# 0.5
mroc.v1 <- roc(as.numeric(datos.vb$RESPONSE), vb.predict.o)
cont <- coords(mroc.v1, .5, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.6
cont <- coords(mroc.v1, .6, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.7
cont <- coords(mroc.v1, .7, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.8
cont <- coords(mroc.v1, .8, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.9
cont <- coords(mroc.v1, .9, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))

################# combined
# 0.5
mroc.v1 <- roc(as.numeric(combined$RESPONSE), combined.predict.o)
cont <- coords(mroc.v1, .5, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.6
cont <- coords(mroc.v1, .6, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.7
cont <- coords(mroc.v1, .7, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.8
cont <- coords(mroc.v1, .8, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))
# 0.9
cont <- coords(mroc.v1, .9, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
sup.tabla.1 <- rbind(sup.tabla.1, sprintf("%.1f", 100*valores))

write.table(sup.tabla.1,file = 'eTable_1', quote = FALSE,sep = "\t")


##################################################
############ Relative Importance #################
##################################################

importance <- datos
modpres <- glm(formula = RESPONSE ~ Onset + SCN1A.score, family = binomial(link = "logit"), data = importance)
dapres<-dominanceAnalysis(modpres)
dominanceMatrix(dapres, type="complete",fit.functions = "r2.m", ordered=TRUE)
plot(dapres, which.graph ="general",fit.function = "r2.m")
pdf(file = "eFigure-1.pdf",width = 4,height = 5) 
plot(dapres, which.graph ="general",fit.function = "r2.m")
dev.off()
dapres


##################################################
################ Random Model ####################
##################################################

aaa <- data.frame(R=datos$RESPONSE, S=datos$SCN1A.score, O=datos$Onset)
bbb <- data.frame(R=datos.va$RESPONSE, S=datos.va$SCN1A.score, O=datos.va$Onset)
ccc <- data.frame(R=datos.vb$RESPONSE, S=datos.vb$SCN1A.score, O=datos.vb$Onset)
random.input <- rbind(aaa,bbb,ccc)
set.seed(123)
new.rows <- sample(nrow(random.input))
random.df <- random.input[new.rows,] 
setenta <- random.df[1:713,]
treinta <- random.df[714:1018,]

rr= datadist(setenta)
options(datadist='rr')
# random traning cohort Modeling
mr<-lrm(R~rcs(O,3) + S, x=TRUE , y=TRUE, setenta)
sr = Score(list('model' = mr), 
           formula=R~1,
           data=setenta,
           plots=c("roc", "calibration"),
           B=1000, trainseeds = 1:1000, summary="IPA")
# random validation
sr.va = Score(list('model' = mr),
             formula=R~1,            
             data=treinta,
             plots=c("roc", "calibration"),
             B=1, summary="IPA")
# Extracting values
sr.data <- as.data.frame(sr$ROC$plotframe)
sr.va.data <- as.data.frame(sr.va$ROC$plotframe)
sr.data$model <- "Random training cohort (70%)"
sr.va.data$model <- "Random validation cohort (30%)"
sr.two <- as.data.frame(rbind(sr.data, sr.va.data))
sr.two$model <- factor(sr.two$model, levels = c("Random training cohort (70%)", "Random validation cohort (30%)"))
# Plot
eFigure2 <- ggplot(data = sr.two, aes(x=(sr.two$FPR), y=sr.two$TPR, colour=sr.two$model),show.legend = FALSE) +
  scale_colour_manual(values=c("dodgerblue3", "cyan")) +
  geom_line(size=0.5) +
  geom_abline() +
  theme_bw() +
  theme(text = element_text(size=8)) +
  labs(x = "1-Specificity",
       y = "Sensitivity",
       colour = element_blank())

ggsave(filename = "eFigure-2.pdf", plot = eFigure2, width = 6.8, height = 5)

random.auc <- rbind(sr$AUC$score, sr.va$AUC$score)
random.auc

# Model performance in Validation Cohorts 1 ans 2 ----------------------------------------------------
m<-lrm(RESPONSE~rcs(Onset,3) + SCN1A.score, x=TRUE , y=TRUE, datos)
s = Score(list('model' = m),
          formula=RESPONSE~1,
          data=datos,
          plots=c("roc", "calibration"),
          B=1000, trainseeds = 1:1000, summary="IPA")

plotCalibration(s, pseudo = FALSE, rug = FALSE, col = "dodgerblue3",)
# Extracting values
modelo = s
auc = modelo$AUC$score$AUC*100
auc.ci1 = modelo$AUC$score$lower*100
auc.ci2 = modelo$AUC$score$upper*100
# IPA 
ipa = modelo$Brier$score$IPA[2]*100
brier = modelo$Brier$score$Brier[1]*100
up = modelo$Brier$contrasts$upper*100*-1
low = modelo$Brier$contrasts$lower*100*-1
ipa.ci1 = up/brier*100
ipa.ci2 = low/brier*100
stats.raw = round(cbind(auc, auc.ci1, auc.ci2, ipa, ipa.ci1, ipa.ci2),2)
stats = cbind("Training", stats.raw)
modelo.results2 = rbind(stats)

# Validation cohort 1
s.va = Score(list('model' = m),
             formula=RESPONSE~1,            
             data=datos.va,
             plots=c("roc", "calibration"),
             B=1, summary="IPA")

# Extracting values
modelo = s.va
auc = modelo$AUC$score$AUC*100
auc.ci1 = modelo$AUC$score$lower*100
auc.ci2 = modelo$AUC$score$upper*100
# IPA 
ipa = modelo$Brier$score$IPA[2]*100
brier = modelo$Brier$score$Brier[1]*100
up = modelo$Brier$contrasts$upper*100*-1
low = modelo$Brier$contrasts$lower*100*-1
ipa.ci1 = up/brier*100
ipa.ci2 = low/brier*100
stats.raw = round(cbind(auc, auc.ci1, auc.ci2, ipa, ipa.ci1, ipa.ci2),2)
stats = cbind("Validation cohort 1", stats.raw)
modelo.results2 = rbind(modelo.results2, stats)

# for vb
s.vb = Score(list('model' = m),
             formula=RESPONSE~1,            
             data=datos.vb,
             plots=c("roc", "calibration"),
             B=1, summary="IPA") 

# Extracting values
modelo = s.vb
auc = modelo$AUC$score$AUC*100
auc.ci1 = modelo$AUC$score$lower*100
auc.ci2 = modelo$AUC$score$upper*100
# IPA 
ipa = modelo$Brier$score$IPA[2]*100
brier = modelo$Brier$score$Brier[1]*100
up = modelo$Brier$contrasts$upper*100*-1
low = modelo$Brier$contrasts$lower*100*-1
ipa.ci1 = up/brier*100
ipa.ci2 = low/brier*100
stats.raw = round(cbind(auc, auc.ci1, auc.ci2, ipa, ipa.ci1, ipa.ci2),2)
stats = cbind("Validation cohort 2", stats.raw)
modelo.results2 = rbind(modelo.results2, stats)
modelo.results2

