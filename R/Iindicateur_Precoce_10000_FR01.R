# Early warning indicators - Stage Risquapi INRA
# Script creation : Fabrice Requier
# Script update & changes : Solene Marion and François Rebaudo
# Last update: 2022-02-14

# --------------------------------------------------------------------------
# Script early warning indicators / 10000 simulations
# --------------------------------------------------------------------------

rm(list=ls())

# 1. Packages
# --------------------------------------------------------------------------
library("lattice")
library("pROC")
library("RColorBrewer")
library("fields")
library("caTools")
library("ggplot2")
library("visreg")

# 2. Dataset
# --------------------------------------------------------------------------
simulRes <- read.table(
  "input_data/resultCombinedStressMiseFormeSoleneTOT10000.csv", 
  sep = ";", dec = ".", header = TRUE
)

simulRes$survival <- as.logical(simulRes$survival)

# simulRes$TotalWeight <- NA
# calcul pour obtenir le poids de la colonie
simulRes$TotalWeight <-
  (
    (simulRes$TotalIHbees + simulRes$TotalForagers)*0.1 + 
    simulRes$TotalDrones*0.2 + 
    simulRes$TotalLarvae*0.0853 + 
    simulRes$TotalDroneLarvae*0.1137333
  )/1000 + 
    simulRes$HoneyEnergyStore 
# ne pas prendre les poids négatifs
simulRes$TotalWeight[simulRes$TotalWeight < 0] <- 0 
# ne pas prendre les poids négatifs
simulRes$HoneyEnergyStore[simulRes$HoneyEnergyStore < 0] <- 0 


# 3. Data analysis
# --------------------------------------------------------------------------


##################################################################################
##################################### Matrice de chaleur #########################
##################################################################################
#############################plusieurs mod??les mais simple#######################################
###### Initialisation de la boucle: pr??paration des objets de sortie
####### param??tres d'analyse:
# ttime_min : jour o?? commencent les analyses
ttime_min <- 100 
# ttime_max : jour o?? finissent les analyses
ttime_max <- 220
# combien de temps plus tard on regarde 
pas <- 10
# Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse
seuil_survival <- 0.99  

####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 485-ttime_max

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- seq(ttime_min,ttime_max,pas)
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

# out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)

#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i <- 1:i_max # seq(1, i_max, 1)
out$nb_survival <- NA
out$nb_survival_n1 <- NA
out$prop_survival_n1 <- NA

out$AUCpop <- NA
out$AUClarve <- NA
out$AUClarvedrone <-NA
out$AUChoneyenergy <- NA
out$AUCtotmite <- NA
out$AUCAFF <- NA
out$AUCLPS <-NA
out$AUCtotevendtoday <- NA
out$AUCDRONE <- NA
out$AUCweight <- NA

out$modpopNRJ <- NA # mod??le selectionner par le stepwise lors de l'interaction population * miel
out$modpopDL <- NA  # mod??le selectionner par le stepwise lors de l'interaction population * drone
out$modLarveDL <- NA  # mod??le selectionner par le stepwise lors de l'interaction larve * drone
out$modLarveNRJ <- NA  # mod??le selectionner par le stepwise lors de l'interaction larve * miel
out$modNRJDL <- NA  # mod??le selectionner par le stepwise lors de l'interaction drone * miel
out$AUCpopNRJ <-NA # population en interation avec miel
out$AUCpopDL <-NA # population en interaction avec drone
out$AUCLarveDL <-NA # larve et drone en interaction
out$AUCLarveNRJ <-NA # larve et miel en interaction
out$AUCNRJDL <-NA # larve et drone en interaction

out$pvaluepop <- NA
out$pvaluelarve <- NA
out$pvaluehoneyenergy <- NA
out$pvaluetotmite <- NA
out$pvalueAFF <- NA
out$pvalueLSP <- NA
out$pvaluetotevendtoday <- NA
out$pvalueweight <- NA

###trouv?? le signe de la pente du glm
out$signepop <- NA
out$signelarve <- NA
out$signelarvedrone <-NA
out$signehoneyenergy <- NA
out$signetotmite <- NA
out$signeAFF<- NA
out$signeLPS <-NA
out$signetotevendtoday <- NA
out$signedrone <- NA
out$signeweight <- NA

### boucle

i<-1

for(i in 1:i_max) {

  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime == time_i & survival == 1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime == time_i + interv_i)[,c("param", "survival")] # donn??e future 
  colnames(futur_i)[2] <- "survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur <- merge(etat_i, futur_i, by = "param") #importation des donn??es survie future dans un m??me tableau
  etat_futur$TotalWeight <-
    (
      (etat_futur$TotalIHbees + etat_futur$TotalForagers)*0.1 + 
      etat_futur$TotalDrones*0.2 + 
      etat_futur$TotalLarvae*0.0853 +
      etat_futur$TotalDroneLarvae*0.1137333/12.78
    )/1000 + etat_futur$HoneyEnergyStore
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"] < seuil_survival){
    
    mod1 <- glm(survival_n1 ~ TotalPopSize, family = binomial, data = etat_futur)
    mod2 <- glm(survival_n1 ~ TotalLarvae,  family = binomial, data = etat_futur)
    mod3 <- glm(survival_n1 ~ TotalDroneLarvae,  family = binomial, data = etat_futur)
    mod4 <- glm(survival_n1 ~ HoneyEnergyStore,  family = binomial, data = etat_futur)
    mod5 <- glm(survival_n1 ~ TotalMites,  family = binomial, data = etat_futur)
    mod6 <- glm(survival_n1 ~ AFF,  family = binomial, data = etat_futur)
    mod7 <- glm(survival_n1 ~ LSP,  family = binomial, data = etat_futur)
    mod8 <- glm(survival_n1 ~ TotalEventsToday,  family = binomial, data = etat_futur)
    mod9 <- glm(survival_n1 ~ TotalDrones, family = binomial, data = etat_futur)
    mod10 <- glm(survival_n1 ~ TotalWeight, family = binomial, data = etat_futur)
    
    mod101 <- glm(survival_n1 ~ TotalPopSize*TotalLarvae, family = binomial, data = etat_futur)
    mod102 <- glm(survival_n1 ~ TotalPopSize*HoneyEnergyStore, family = binomial, data = etat_futur)
    mod103 <- glm(survival_n1 ~ TotalPopSize*TotalDroneLarvae, family = binomial, data = etat_futur)
    mod104 <- glm(survival_n1 ~ TotalLarvae*HoneyEnergyStore, family = binomial, data = etat_futur)
    mod105 <- glm(survival_n1 ~ TotalLarvae*TotalDroneLarvae, family = binomial, data = etat_futur)
    mod106 <- glm(survival_n1 ~ HoneyEnergyStore*TotalDroneLarvae, family = binomial, data = etat_futur)
    #mod??le selectionner par le stepwise
    modele_selectionne1 <- step(mod101)
    modele_selectionne2 <- step(mod102)
    modele_selectionne3 <- step(mod103)
    modele_selectionne4 <- step(mod104)
    modele_selectionne5 <- step(mod105)
    modele_selectionne6 <- step(mod106)
    modeleselec1 <-summary(modele_selectionne1)
    modeleselec2 <-summary(modele_selectionne2)
    modeleselec3 <-summary(modele_selectionne3)
    modeleselec4 <-summary(modele_selectionne4)
    modeleselec5 <-summary(modele_selectionne5)
    modeleselec6 <-summary(modele_selectionne6)
    modeleselec21<-as.character(modeleselec1$call)
    modeleselec22<-as.character(modeleselec2$call)
    modeleselec23<-as.character(modeleselec3$call)
    modeleselec24<-as.character(modeleselec4$call)
    modeleselec25<-as.character(modeleselec5$call)
    modeleselec26<-as.character(modeleselec6$call)
    out[i, "modint"] <- modeleselec21 [2]
    out[i, "modpopNRJ"] <- modeleselec22 [2]
    out[i, "modpopDL"] <- modeleselec23 [2]
    out[i, "modLarveDL"] <- modeleselec24 [2]
    out[i, "modLarveNRJ"] <- modeleselec25 [2]
    out[i, "modNRJDL"] <- modeleselec26 [2]
    
    roc1 <- roc(etat_futur[names(fitted(mod1)), "survival_n1"] ~ fitted(mod1))
    roc2 <- roc(etat_futur[names(fitted(mod2)), "survival_n1"] ~ fitted(mod2))
    roc3 <- roc(etat_futur[names(fitted(mod3)), "survival_n1"] ~ fitted(mod3))
    roc4 <- roc(etat_futur[names(fitted(mod4)), "survival_n1"] ~ fitted(mod4))
    roc5 <- roc(etat_futur[names(fitted(mod5)), "survival_n1"] ~ fitted(mod5))
    roc6 <- roc(etat_futur[names(fitted(mod6)), "survival_n1"] ~ fitted(mod6))
    roc7 <- roc(etat_futur[names(fitted(mod7)), "survival_n1"] ~ fitted(mod7))
    roc8 <- roc(etat_futur[names(fitted(mod8)), "survival_n1"] ~ fitted(mod8))
    roc9 <- roc(etat_futur[names(fitted(mod9)), "survival_n1"] ~ fitted(mod9))
    roc10 <- roc(etat_futur[names(fitted(mod10)), "survival_n1"] ~ fitted(mod10))
    roc101 <- roc(etat_futur[names(fitted(mod101)), "survival_n1"] ~ fitted(mod101)) #MS = Modele Selectionn??
    roc102 <- roc(etat_futur[names(fitted(mod102)), "survival_n1"] ~ fitted(mod102)) #MS = Modele Selectionn??
    roc103 <- roc(etat_futur[names(fitted(mod103)), "survival_n1"] ~ fitted(mod103)) #MS = Modele Selectionn??
    roc104 <- roc(etat_futur[names(fitted(mod104)), "survival_n1"] ~ fitted(mod104)) #MS = Modele Selectionn??
    roc105 <- roc(etat_futur[names(fitted(mod105)), "survival_n1"] ~ fitted(mod105)) #MS = Modele Selectionn??
    roc106 <- roc(etat_futur[names(fitted(mod106)), "survival_n1"] ~ fitted(mod106)) #MS = Modele Selectionn??
    
    AUC1 <- auc(roc1)
    AUC2 <- auc(roc2)
    AUC3 <- auc(roc3)
    AUC4 <- auc(roc4)
    AUC5 <- auc(roc5)
    AUC6 <- auc(roc6)
    AUC7 <- auc(roc7)
    AUC8 <- auc(roc8)
    AUC9 <- auc(roc9)
    AUC10 <- auc(roc10)
    AUCMS1 <- auc(roc101)
    AUCMS2 <- auc(roc102)
    AUCMS3 <- auc(roc103)
    AUCMS4 <- auc(roc104)
    AUCMS5 <- auc(roc105)
    AUCMS6 <- auc(roc106)
    
    out[i,"AUCpop"] <- AUC1
    out[i,"AUClarve"] <- AUC2
    out[i,"AUClarvedrone"] <- AUC3
    out[i,"AUChoneyenergy"] <- AUC4
    out[i,"AUCtotmite"] <- AUC5
    out[i,"AUCAFF"] <- AUC6
    out[i,"AUCLPS"] <- AUC7
    out[i,"AUCtotevendtoday"] <- AUC8
    out[i, "AUCDRONE"] <- AUC9
    out[i,"AUCweight"] <- AUC10
    
    out[i,"AUCpopNRJ"] <- AUCMS2
    out[i,"AUCpopDL"] <- AUCMS3
    out[i,"AUCLarveNRJ"] <- AUCMS4
    out[i,"AUCLarveDL"] <- AUCMS5
    out[i,"AUCNRJDL"] <- AUCMS6
    
    out[i,"pvaluepop"] <- coef(summary(mod1))[2,4]
    out[i,"pvaluelarve"] <- coef(summary(mod2))[2,4]
    out[i,"pvaluehoneyenergy"] <- coef(summary(mod4))[2,4]
    out[i,"pvaluetotmite"] <- coef(summary(mod5))[2,4]
    out[i,"pvalueAFF"] <- coef(summary(mod6))[2,4]
    out[i,"pvalueLSP"] <- coef(summary(mod7))[2,4]
    out[i,"pvaluetotevendtoday"] <- coef(summary(mod8))[2,4]
    out[i,"pvalueweight"] <- coef(summary(mod10))[2,4]
    
    
    out[i,"signepop"] <- (sign(coefficients(mod1)[2]))
    out[i,"signelarve"] <- (sign(coefficients(mod2)[2]))
    out[i,"signelarvedrone"] <- (sign(coefficients(mod3)[2]))
    out[i,"signehoneyenergy"] <- (sign(coefficients(mod4)[2]))
    out[i,"signetotmite"] <- (sign(coefficients(mod5)[2]))
    out[i,"signeAFF"] <- (sign(coefficients(mod6)[2]))
    out[i,"signeLPS"] <- (sign(coefficients(mod7)[2]))
    out[i,"signetotevendtoday"] <- (sign(coefficients(mod8)[2]))
    out[i,"signedrone"] <- (sign(coefficients(mod9)[2]))
    out[i,"signeweight"] <- (sign(coefficients(mod10)[2]))
  }
} 

save(out, file = "./input_data/out.RData")


###Faire les matrices de chaleurs en se basant sur l'AUC

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256) # la base de couleur utilis??e pour la matrice de chaleur
par(mfrow = c(1,1))
#proportion survie
outmatprop_survival_n1 <- matrix(out$prop_survival_n1, length(interv), length(ttime_n), byrow = TRUE)
fields::image.plot(interv, ttime_n, outmatprop_survival_n1, col=rev(cols), main = "Proportion de survie", zlim = c(0.5, 1))
contour(interv, ttime_n, outmatprop_survival_n1, levels = seq(0, 1, by = 0.05), add = TRUE, col = "black")


#matrice AUC

#pop
outmatpop <- matrix(out$AUCpop, length(interv), length(ttime_n), byrow = TRUE)
fields::image.plot(interv, ttime_n, outmatpop, col = rev(cols), main ="AUC taille population",zlim = c(0.5, 1))
contour(interv,ttime_n, outmatpop, levels = seq(0.5, 1, by = 0.05), add = TRUE, col = "black")
#AUClarve
outmatlarve<- matrix(out$AUClarve,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatlarve, col=rev(cols), main ="Couvain ouvri??re",zlim=c(0.5, 1))
contour(interv,ttime_n, outmatlarve, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUC Honey Energy Store
outmatpvaluehoneyenergy<- matrix(out$AUChoneyenergy,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatpvaluehoneyenergy, col=rev(cols), main ="AUC Honey Energy Store", zlim=c(0.5, 1))
contour(interv,ttime_n, outmatpvaluehoneyenergy, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUC larve drones
outmatDlarve<- matrix(out$AUClarvedrone,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatDlarve, col=rev(cols), main ="AUC Larve Drone",zlim=c(0.5, 1) )
contour(interv,ttime_n, outmatDlarve, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUCdrone
outmatdrone<- matrix(out$AUCDRONE,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatdrone, col=rev(cols), main ="AUC M??le adulte", zlim=c(0.5, 1))
contour(interv,ttime_n, outmatdrone, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUC mite
outmatmite<- matrix(out$AUCtotmite,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatmite, col=rev(cols), main ="AUC mite", zlim=c(0.5, 1))
contour(interv,ttime_n, outmatmite, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUC AFF
outmatAFF<- matrix(out$AUCAFF,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatAFF, col=rev(cols), main ="AUC AFF", zlim=c(0.5, 1))
contour(interv,ttime_n, outmatAFF, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUCLPS
outmatLPS<- matrix(out$AUCLPS,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatLPS, col=rev(cols), main ="AUC LSP",zlim=c(0.5, 1))
contour(interv,ttime_n, outmatLPS, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUCtotevendtoday
outmatevend<- matrix(out$AUCtotevendtoday,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatevend, col=rev(cols), main ="AUC Total evend day", zlim=c(0.5, 1))
contour(interv,ttime_n, outmatevend, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUCweight
outmatweight<- matrix(out$AUCweight,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatweight, col=rev(cols), main ="AUC Poids de la colonie",xlab="Pr??diction (en jours)", ylab="Jours d'observation", zlim=c(0.4, 1))
contour(interv,ttime_n, outmatweight, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")

#matrice pop*NRJ
matAUCpopNRJ<- matrix(out$AUCpopNRJ,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matAUCpopNRJ, col=rev(cols), main ="AUC interaction pop NRJ", zlim = c(0.5, 1))
contour(interv,ttime_n, matAUCpopNRJ, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#matrice pop*DL
matAUCpopDL<- matrix(out$AUCpopDL,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matAUCpopDL, col=rev(cols), main ="AUC interaction pop DL", zlim = c(0.5, 1))
contour(interv,ttime_n, matAUCpopDL, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#matrice larve*DL
matAUCLarveDL<- matrix(out$AUCLarveDL,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matAUCLarveDL, col=rev(cols), main ="AUC interaction Larve DL", zlim = c(0.5, 1))
contour(interv,ttime_n, matAUCLarveDL, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#matrice larve*NRJ
matAUCLarveNRJ<- matrix(out$AUCLarveNRJ,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matAUCLarveNRJ, col=rev(cols), main ="AUC interaction Larve NRJ", zlim = c(0.5, 1))
contour(interv,ttime_n, matAUCLarveNRJ, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#matrice DL*NRJ
matAUCNRJDL<- matrix(out$AUCNRJDL,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matAUCNRJDL, col=rev(cols), main ="AUC interaction DL NRJ", zlim = c(0.5, 1))
contour(interv,ttime_n, matAUCNRJDL, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")




###matrice p value

#pop 
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
outmat<- matrix(out$AUCpop,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmat, col=rev(cols), main ="")
contour(interv,ttime_n, outmat, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUClarve
outmatpvaluelarve<- matrix(out$pvaluelarve,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatpvaluelarve, col=rev(cols), main ="p value larve")
contour(interv,ttime_n, outmatpvaluelarve, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUC honey
outmathoney<- matrix(out$AUChoneyenergy,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmathoney, col=rev(cols), main ="p value honey")
contour(interv,ttime_n, outmathoney, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUC mite
outmatpvaluetotmite<- matrix(out$pvaluetotmite,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatpvaluetotmite, col=rev(cols), main ="p value mite")
contour(interv,ttime_n, outmatpvaluetotmite, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUC AFF
outmatpvalueAFF<- matrix(out$pvalueAFF,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatpvalueAFF, col=rev(cols), main ="p value aff")
contour(interv,ttime_n, outmatpvalueAFF, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUCLPS
outmatpvalueLSP<- matrix(out$pvalueLSP,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatpvalueLSP, col=rev(cols), main ="p value lsp")
contour(interv,ttime_n, outmatpvalueLSP, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUCtotevendtoday
outmapvaluetotevendtoday<- matrix(out$pvaluetotevendtoday,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmapvaluetotevendtoday, col=rev(cols), main ="p value evend")
contour(interv,ttime_n, outmapvaluetotevendtoday, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUCweight
outmatpvalueweight<- matrix(out$pvalueweight,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatpvalueweight, col=rev(cols), main ="p value Poids")
contour(interv,ttime_n, outmatpvalueweight, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")

###matrice signe

#signe larve
matlarvesigne<- matrix(out$signelarve,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matlarvesigne, col=rev(cols), main ="Larve Signe")
contour(interv,ttime_n, matlarvesigne, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#signe pop
matpopsigne<- matrix(out$signepop,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matpopsigne, col=rev(cols), main ="pop signe")
contour(interv,ttime_n, matpopsigne, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#signe Honey Energy Store
mathoneysigne<- matrix(out$signehoneyenergy,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, mathoneysigne, col=rev(cols), main ="Honey Energy Store signe")
contour(interv,ttime_n, mathoneysigne, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#Signe larve drones
matlarvedronesigne<- matrix(out$signelarvedrone,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matlarvedronesigne, col=rev(cols), main ="Larve Drone signe")
contour(interv,ttime_n, matlarvedronesigne, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#signe mite
matmitesigne<- matrix(out$signetotmite,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matmitesigne, col=rev(cols), main ="mite signe")
contour(interv,ttime_n, matmitesigne, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#signe AFF
matAFFsigne<- matrix(out$signeAFF,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matAFFsigne, col=rev(cols), main ="AFF signe")
contour(interv,ttime_n, matAFFsigne, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#signe LPS
matLPSsigne<- matrix(out$signeLPS,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matLPSsigne, col=rev(cols), main ="LSP signe")
contour(interv,ttime_n, matLPSsigne, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#signe totevendtoday
matevendsigne<- matrix(out$signetotevendtoday,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matevendsigne, col=rev(cols), main ="Total evend day signe")
contour(interv,ttime_n, matevendsigne, levels = seq(-1, 1, by = 0.5),add = TRUE, col = "black")
#signe drone
matdronesigne<- matrix(out$signedrone,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matdronesigne, col=rev(cols), main ="Total evend day signe")
contour(interv,ttime_n, matdronesigne, levels = seq(-1, 1, by = 0.5),add = TRUE, col = "black")
#signe Poids
matweightsigne <- matrix(out$signeweight,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matweightsigne, col=rev(cols), main ="Poids signe")
contour(interv,ttime_n, matweightsigne, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")

