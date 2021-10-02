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
datos <- as.data.frame(read.table("SCN1A_Training_Cohort.txt", header = TRUE, sep = "\t"))
#datos <- as.data.frame(read.table("No_GEFS+_related.txt", header = TRUE, sep = "\t"))
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
################################################## -------------------------

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
ggsave(filename = "Figure-2.pdf", plot = Figure2, width = 5, height = 6)


##################################################
################# Main modeling ##################
##################################################

# SCN1A score & Onset model (Final) -------------------------------------------------------------
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

# SCN1A Genetic score Model --------------------------------------------------------
mp<-lrm(r~p, x=TRUE , y=TRUE)
sp = Score(list('model' = mp), 
           formula=RESPONSE~1,
           data=datos,
           plots=c("roc", "calibration"),
           B=1000, trainseeds = 1:1000, summary="IPA")

# Extracting value
modelo = sp
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
stats = cbind("SCN1A score", stats.raw)
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
################################################## ---------------------------------------------------

plotCalibration(s, pseudo = FALSE, rug = FALSE, col = "dodgerblue3")
plotCalibration(sp, pseudo = FALSE, rug = FALSE, col = "gold")
plotCalibration(so, pseudo = FALSE, rug = FALSE, col = "darkorange1")
plotCalibration(sc, pseudo = FALSE, rug = FALSE, col = "forestgreen")
plotCalibration(sv, pseudo = FALSE, rug = FALSE, col = "purple")

# ggplot version
s.data <- as.data.frame(s$ROC$plotframe)
so.data <- as.data.frame(so$ROC$plotframe)
sc.data <- as.data.frame(sc$ROC$plotframe)
sv.data <- as.data.frame(sv$ROC$plotframe)
sp.data <- as.data.frame(sp$ROC$plotframe)
cero <- data.frame(model ="SCN1A genetic score", risk=1,TPR=0,FPR=0)
sp.data <- rbind(sp.data,cero)
s.data$model <- "SCN1A genetic score + Onset"
so.data$model <- "Onset-only"
sc.data$model <- "CADD + Onset"
sv.data$model <- "REVEL + Onset"
sp.data$model <- "SCN1A genetic score"
# TPR is sensitivity(y-axis), 1-FPR is =1-Specificity (x-axis)
s.all <- as.data.frame(rbind(s.data, so.data, sc.data, sv.data, sp.data))
s.all$model <- factor(s.all$model, levels = c("SCN1A genetic score + Onset", "Onset-only","CADD + Onset","REVEL + Onset", "SCN1A genetic score"))
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

ggsave(filename = "Figure-3.pdf", plot = Figure3, width = 5, height = 5)

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

#Combined onset.only
combined.onset <- data.frame(o=combined$Onset)
combined.predict.o <- predict(mo, combined.onset, type="fitted.ind")

##################################################
#################### FIGURE 5 ####################
##################################################

# Validation Cohort Plot BARPLOT AUSTRALIA ---------------------------------------------
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

# Validation Cohort Plot BARPLOT BELGIUM ---------------------------------------------
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

# Validation Cohort Density Plot AUSTRALIA ------------------------------------------
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

# Validation Cohort Density Plot BELGIUM------------------------------------------
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

#########################################################################
############ Supplementary_Material_Cohorts_composition #################
#########################################################################

# Version 2 - by Cohort and Phenotype
Training <- as.data.frame(datos$Phenotype)
names(Training)[names(Training) == "datos$Phenotype"] <- "Phenotype"
Training$set <- "Training cohort"
Validation1 <- as.data.frame(datos.va$Phenotype)
names(Validation1)[names(Validation1) == "datos.va$Phenotype"] <- "Phenotype"
Validation1$set <- "Validation cohort 1"

Validation2 <- as.data.frame(datos.vb$Phenotype)
names(Validation2)[names(Validation2) == "datos.vb$Phenotype"] <- "Phenotype"
Validation2$set <- "Validation cohort 2"

s1 <- ggplot(data = Training, aes(x=Phenotype, fill = Phenotype, y = ((..count..)/sum(..count..)*100))) +
  scale_fill_manual(labels=c("Dravet Syndrome", "GEFS+"),
                    breaks=c("Dravet","GEFS+"),
                    values=c("mediumpurple4","grey48")) +
  geom_bar(show.legend = FALSE) +
  ylim(0, 100) +
  theme_bw() +
  theme(text = element_text(size=8)) +
  labs(title="Training cohort (n=745)",
       x = element_blank(),
       y = "Patients (%)", fill = element_blank())
s2 <- ggplot(data = Validation1, aes(x=Phenotype, fill = Phenotype, y = ((..count..)/sum(..count..)*100))) +
  scale_fill_manual(labels=c("Dravet Syndrome", "GEFS+"),
                    breaks=c("Dravet","GEFS+"),
                    values=c("mediumpurple4","grey48")) +
  geom_bar(show.legend = FALSE) +
  ylim(0, 100) +
  theme_bw() +
  theme(text = element_text(size=8)) +
  labs(title="Validation cohort 1 (n=203)",
       x = element_blank(),
       y = "Patients (%)", fill = element_blank())
s3 <- ggplot(data = Validation2, aes(x=Phenotype, fill = Phenotype, y = ((..count..)/sum(..count..)*100))) +
  scale_fill_manual(labels=c("Dravet Syndrome", "GEFS+"),
                    breaks=c("Dravet","GEFS+"),
                    values=c("mediumpurple4","grey48")) +
  geom_bar(show.legend = FALSE) +
  ylim(0, 100) +
  theme_bw() +
  theme(text = element_text(size=8)) +
  labs(title="Validation cohort 2 (n=72)",
       x = element_blank(),
       y = "Patients (%)", fill = element_blank())

lay <- rbind(c(1,2,3))
Cohort.composition <- grid.arrange(grobs = list(s1,s2,s3), layout_matrix = lay)
ggsave(filename = "Supplementary_Material_Cohort_composition.pdf", plot = Cohort.composition, width = 6.6, height = 3.3)

#########################################################################
############ Supplementary_Material_Relative_Importance #################
#########################################################################

importance <- datos
modpres <- glm(formula = RESPONSE ~ Onset + SCN1A.score, family = binomial(link = "logit"), data = importance)
dapres<-dominanceAnalysis(modpres)
dominanceMatrix(dapres, type="complete",fit.functions = "r2.m", ordered=TRUE)
rip<-plot(dapres, which.graph ="general",fit.function = "r2.m")
dapres
ggsave(filename = "Supplementary_Material_Relative_Importance", plot = rip, width = 5, height = 4)

####################################################################
#################### Supplementary_Material_Table_1 ##################
####################################################################

# Comparison between SCN1A genetic score + Onset and Onset only models across different thresholds. 

#traning.c
traning.c <- datos
traning.onset.pg <- data.frame(o=traning.c$Onset,p=traning.c$SCN1A.score)
traning.onset <- data.frame(o=traning.c$Onset)
traning.pg <- data.frame(p=traning.c$SCN1A.score)
traning.predict <- predict(m, traning.onset.pg, type="fitted.ind")
traning.predict.o <- predict(mo, traning.onset, type="fitted.ind")
traning.predict.p <- predict(mp, traning.pg, type="fitted.ind")
tabla.supp <- c("Sensitivity","ci1","ci2","Specificity","ci1","ci2","PPV","ci1","ci2","NPV","ci1","ci2","Accuracy","ci1","ci2")

mroc.v1 <- roc(as.numeric(traning.c$RESPONSE), traning.predict)

cont <- coords(mroc.v1, .5, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .6, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .7, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .8, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .9, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

# onset.only
mroc.v1 <- roc(as.numeric(traning.c$RESPONSE), traning.predict.o)
cont <- coords(mroc.v1, .5, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .6, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .7, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .8, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .9, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

# SCN1A genetic score
mroc.v1 <- roc(as.numeric(traning.c$RESPONSE), traning.predict.p)
cont <- coords(mroc.v1, .5, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .6, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .7, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .8, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

cont <- coords(mroc.v1, .9, "threshold", ret=c("tp","fp","fn","tn"))
tabla.c <- as.table(matrix(c(cont$tp,cont$fp,cont$fn,cont$tn), nrow = 2, byrow = TRUE))
rval <- epi.tests(tabla.c, conf.level = 0.95)
valores <- cbind(rval$detail$se, rval$detail$sp, rval$detail$pv.pos,rval$detail$pv.neg, rval$detail$diag.ac)
tabla.supp <- rbind(tabla.supp, sprintf("%.1f", 100*valores))

write.table(tabla.supp,file = "Supplementary_Material_Table_1",quote = FALSE, sep = "\t")

# END

