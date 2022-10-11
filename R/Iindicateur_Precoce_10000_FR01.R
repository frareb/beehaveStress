# Early warning indicators - Stage Risquapi INRA
# Script creation : Fabrice Requier
# Script update & changes : Solene Marion and François Rebaudo
# Last update: 2022-02-14

# --------------------------------------------------------------------------
# Script early warning indicators / 10000 simulations
# --------------------------------------------------------------------------

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
  "resultCombinedStressMiseFormeSoleneTOT10000.csv", 
  sep = ";", dec = ".", header = TRUE
)
str(simulRes)
simulRes$survival <- as.logical(simulRes$survival)
unique(simulRes$DateCollapse)

# 3. Data analysis
# --------------------------------------------------------------------------
# Percentage of alive and dead colonies
# dead colonies
vecDead <- unique(simulRes[simulRes$survival == FALSE, "param"])
deadColoniesP <- length(vecDead)/length(unique(simulRes$param))
# alive colonies
vecAlive = unique(simulRes$param)[is.na(match(unique(simulRes$param),vecDead))]
aliveColoniesP <- length(vecAlive)/length(unique(simulRes$param))
# dead / alive colonies percentage
print(data.frame(deadColoniesP , aliveColoniesP))





##### Courbe de survie des colonies en fonction du temps 
#obtention du tableau de survie en fonction du temps
simulRes$survie=ifelse(simulRes$survival=="1", "survive", "meurt")
tableau <- table (simulRes$ttime, simulRes$survie)
#charg?? le tableau d??j?? fait
tableau <- read.table(
  "/Users/solene/Desktop/out graph ect/time_meurt_vie.csv", 
  sep = ";", dec = ".", header = TRUE
)
plot(
  x = tableau$Time, 
  y = tableau$survive, 
  type = "l", 
  main = "Survie des colonies en fonction du temps", 
  xlab = "Temps (en jours)", ylab = "Nombre de colonie en vie", lwd = 5
)

#####Graphique, ??volution du nombre de colonies morte en fonction du temps

plot(simulRes$ttime, main = "nombre de colonie morte en fonction du temps", xlab ="temps", ylab ="nombre de colonie morte")

out <- read.table("/Users/solene/Desktop/out graph ect/fichier excel out/SUPER OUT COMPLET.csv", sep=";",dec=".",h=T)


##################################################################################
##################################### Matrice de chaleur #########################
##################################################################################

#############################1 seul mod??le #######################################

####### param??tres d'analyse:
# ttime_min : jour o?? commencent les analyses
ttime_min<-100 
#ttime_max : jour o?? finissent les analyses
ttime_max<-300 
# combien de temps plus tard on regarde 
pas <- 10

seuil_survival <- 0.95  #Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse
seuil_AUCsigne <- 0.80
####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 485-ttime_max

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- seq(ttime_min,ttime_max,pas)
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
head(out)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i<-seq(1,i_max,1)
out$nb_survival <- NA
out$nb_survival_n1 <- NA
out$prop_survival_n1 <- NA
#out$AICout <- NA
out$AUCpop <- NA
out$pvalue <- NA

### boucle

#i<-4
i <-1
for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod1 <- glm(survival_n1 ~ TotalPopSize, family = binomial, data = etat_futur)
    #out[i,"AICout"] <- AIC(mod1)
    roc1 <- roc(etat_futur[names(fitted(mod1)), "survival_n1"] ~ fitted(mod1))
    AUC1 <- auc(roc1)
    out[i,"AUCpop"] <- AUC1
    out[i,"pvalue"] <- coef(summary(mod1))[2,4] # extraction de la p value
  }
}

#matrice AUC
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
outmat<- matrix(out$AUCpop,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmat, rev(col=cols), main ="AUC")
contour(interv,ttime_n, outmat, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")

#############################plusieurs mod??les mais simple#######################################


####### param??tres d'analyse:
# ttime_min : jour o?? commencent les analyses
ttime_min<-100 
#ttime_max : jour o?? finissent les analyses
ttime_max<-220
# combien de temps plus tard on regarde 
pas <- 10

seuil_survival <- 0.99  #Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse

####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 485-ttime_max

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- seq(ttime_min,ttime_max,pas)
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
head(out)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i<-seq(1,i_max,1)
out$nb_survival <- NA
out$nb_survival_n1 <- NA
out$prop_survival_n1 <- NA
out$AUCpop <- NA
out$AUClarve <- NA
out$AUClarvedrone <-NA
out$AUChoneyenergy <- NA
out$AUCtotmite <- NA
out$AUCAFF<- NA
out$AUCLPS <-NA
out$AUCtotevendtoday <- NA
out$AUCDRONE <- NA
#out$pvaluepop <- NA
#out$pvaluelarvedrone <- NA
#out$pvaluehoneyenergy <- NA
#out$pvaluetotmite <- NA
#out$pvalueAFF <- NA
#out$pvalueLSP <- NA
#out$pvaluetotevendtoday <- NA

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



### boucle

i<-1

for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod1 <- glm(survival_n1 ~ TotalPopSize, family = binomial, data = etat_futur)
    mod2 <- glm(survival_n1 ~ TotalLarvae,  family = binomial, data = etat_futur)
    mod3 <- glm(survival_n1 ~ TotalDroneLarvae,  family = binomial, data = etat_futur)
    mod4 <- glm(survival_n1 ~ HoneyEnergyStore,  family = binomial, data = etat_futur)
    mod5 <- glm(survival_n1 ~ TotalMites,  family = binomial, data = etat_futur)
    mod6 <- glm(survival_n1 ~ AFF,  family = binomial, data = etat_futur)
    mod7 <- glm(survival_n1 ~ LSP,  family = binomial, data = etat_futur)
    mod8 <- glm(survival_n1 ~ TotalEventsToday,  family = binomial, data = etat_futur)
    mod9 <- glm(survival_n1 ~ TotalDrones, family = binomial, data = etat_futur)
    roc1 <- roc(etat_futur[names(fitted(mod1)), "survival_n1"] ~ fitted(mod1))
    roc2 <- roc(etat_futur[names(fitted(mod2)), "survival_n1"] ~ fitted(mod2))
    roc3 <- roc(etat_futur[names(fitted(mod3)), "survival_n1"] ~ fitted(mod3))
    roc4 <- roc(etat_futur[names(fitted(mod4)), "survival_n1"] ~ fitted(mod4))
    roc5 <- roc(etat_futur[names(fitted(mod5)), "survival_n1"] ~ fitted(mod5))
    roc6 <- roc(etat_futur[names(fitted(mod6)), "survival_n1"] ~ fitted(mod6))
    roc7 <- roc(etat_futur[names(fitted(mod7)), "survival_n1"] ~ fitted(mod7))
    roc8 <- roc(etat_futur[names(fitted(mod8)), "survival_n1"] ~ fitted(mod8))
    roc9 <- roc(etat_futur[names(fitted(mod9)), "survival_n1"] ~ fitted(mod9))
    AUC1 <- auc(roc1)
    AUC2 <- auc(roc2)
    AUC3 <- auc(roc3)
    AUC4 <- auc(roc4)
    AUC5 <- auc(roc5)
    AUC6 <- auc(roc6)
    AUC7 <- auc(roc7)
    AUC8 <- auc(roc8)
    AUC9 <- auc(roc9)
    out[i,"AUCpop"] <- AUC1
    out[i,"AUClarve"] <- AUC2
    out[i,"AUClarvedrone"] <- AUC3
    out[i,"AUChoneyenergy"] <- AUC4
    out[i,"AUCtotmite"] <- AUC5
    out[i,"AUCAFF"] <- AUC6
    out[i,"AUCLPS"] <- AUC7
    out[i,"AUCtotevendtoday"] <- AUC8
    out[i, "AUCDRONE"] <- AUC9
    #out[i,"pvaluepop"] <- coef(summary(mod1))[2,4]
    #out[i,"pvaluelarvae"] <- coef(summary(mod2))[2,4]
    #out[i,"pvaluehoneyenergy"] <- coef(summary(mod4))[2,4]
    #out[i,"pvaluetotmite"] <- coef(summary(mod5))[2,4]
    #out[i,"pvalueAFF"] <- coef(summary(mod6))[2,4]
    #out[i,"pvalueLSP"] <- coef(summary(mod7))[2,4]
    #out[i,"pvaluetotevendtoday"] <- coef(summary(mod8))[2,4]
    out[i,"signepop"] <- (sign(coefficients(mod1)[2]))
    out[i,"signelarve"] <- (sign(coefficients(mod2)[2]))
    out[i,"signelarvedrone"] <- (sign(coefficients(mod3)[2]))
    out[i,"signehoneyenergy"] <- (sign(coefficients(mod4)[2]))
    out[i,"signetotmite"] <- (sign(coefficients(mod5)[2]))
    out[i,"signeAFF"] <- (sign(coefficients(mod6)[2]))
    out[i,"signeLPS"] <- (sign(coefficients(mod7)[2]))
    out[i,"signetotevendtoday"] <- (sign(coefficients(mod8)[2]))
    out[i,"signedrone"] <- (sign(coefficients(mod9)[2]))
  }
} 


###Faire les matrices de chaleurs en se basant sur l'AUC
#matrice AUC
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256) # la base de couleur utilis??e pour la matrice de chaleur
par(mfrow = c(1,1))
#proportion survie
outmatprop_survival_n1<- matrix(out$prop_survival_n1,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatprop_survival_n1, col=rev(cols), main = "Proprotion de survie", zlim=c(0.5, 1))
contour(interv,ttime_n, outmatprop_survival_n1, levels = seq(0, 1, by = 0.05),add = TRUE, col = "black")
par(mfrow = c(3,3))
#AUClarve
outmatlarve<- matrix(out$AUClarve,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatlarve, col=rev(cols), main ="Couvain ouvri??re",zlim=c(0.5, 1))
contour(interv,ttime_n, outmatlarve, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#pop
outmatpop<- matrix(out$AUCpop,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatpop, col=rev(cols), main ="AUC taille population",zlim=c(0.5, 1))
contour(interv,ttime_n, outmatpop, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
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
contour(interv,ttime_n, outmatdrone, contourr,add = TRUE, col = "black")
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

###matrice p value

#pop 
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
outmat<- matrix(out$AUCpop,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmat, col=rev(cols), main ="")
contour(interv,ttime_n, outmat, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUClarve
outmatpvaluelarvae<- matrix(out$pvaluelarvae,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatpvaluelarvae, col=rev(cols), main ="p value larve")
contour(interv,ttime_n, outmatpvaluelarvae, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
#AUC larve drones
outmatDlarve<- matrix(out$AUClarvedrone,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, outmatDlarve, col=rev(cols), main ="p value larve drone")
contour(interv,ttime_n, outmatDlarve, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
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

###matrice signe
par(mfrow = c(3,3))
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

####poids
simulRes$TotalWeight <- NA

simulRes$TotalWeight <-((simulRes$TotalIHbees+simulRes$TotalForagers)*0.1 + simulRes$TotalDrones*0.2+simulRes$TotalLarvae*0.0853+simulRes$TotalDroneLarvae*0.1137333)/1000 +simulRes$HoneyEnergyStore #calcul pour obtenir le poids de la colonie

plot(simulRes$ttime, simulRes$TotalWeight)
simulRes$TotalWeight[simulRes$TotalWeight<0] = 0 #ne pas prendre les poids n??gatifs
simulRes$HoneyEnergyStore[simulRes$HoneyEnergyStore<0]=0 #ne pas prendre les poids n??gatifs

####### param??tres d'analyse:
# ttime_min : jour o?? commencent les analyses
ttime_min<-100 
#ttime_max : jour o?? finissent les analyses
ttime_max<-220 
# combien de temps plus tard on regarde 
pas <- 10

seuil_survival <- 0.99  #Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse

####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 485-ttime_max

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- seq(ttime_min,ttime_max,pas)
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
head(out)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i<-seq(1,i_max,1)
out$nb_survival <- NA
out$nb_survival_n1 <- NA
out$prop_survival_n1 <- NA

out$AUCweight <- NA
out$signeweight <- NA

### boucle

#i<-4
i <-1
for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  etat_futur$TotalWeight <-((etat_futur$TotalIHbees+etat_futur$TotalForagers)*0.1 + etat_futur$TotalDrones*0.2+etat_futur$TotalLarvae*0.0853+
                              etat_futur$TotalDroneLarvae*0.1137333/12.78)/1000 +etat_futur$HoneyEnergyStore
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod30 <- glm(survival_n1 ~ TotalWeight, family = binomial, data = etat_futur)
    roc30 <- roc(etat_futur[names(fitted(mod30)), "survival_n1"] ~ fitted(mod30))
    AUC30 <- auc(roc30)
    out[i,"AUCweight"] <- AUC30
    out[i,"signeweight"] <- (sign(coefficients(mod30)[2]))
  }
}

outmatweight<- matrix(out$AUCweight,length(interv),length(ttime_n),byrow=T)

fields::image.plot(interv, ttime_n, outmatweight, col=rev(cols), main ="Poids de la colonie",xlab="Pr??diction (en jours)", ylab="Jours d'observation", zlim=c(0.4, 1))
contour(interv,ttime_n, outmatweight, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")





################################## Mod??le en interaction ##################################
###############interaction double

####### param??tres d'analyse:
# ttime_min : jour o?? commencent les analyses
ttime_min<-100 
#ttime_max : jour o?? finissent les analyses
ttime_max<-220
# combien de temps plus tard on regarde 
pas <- 10

seuil_survival <- 0.99  #Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse

####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 485-ttime_max

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- seq(ttime_min,ttime_max,pas)
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i<-seq(1,i_max,1)
out$nb_survival <- NA
out$nb_survival_n1 <- NA
out$prop_survival_n1 <- NA
out$modpopNRJ <- NA # mod??le selectionner par le stepwise lors de l'interaction population * miel
out$modpopDL <- NA  # mod??le selectionner par le stepwise lors de l'interaction population * drone
out$modLarveDL <- NA  # mod??le selectionner par le stepwise lors de l'interaction larve * drone
out$modLarveNRJ <- NA  # mod??le selectionner par le stepwise lors de l'interaction larve * miel
out$modNRJDL <- NA  # mod??le selectionner par le stepwise lors de l'interaction drone * miel
out$AUCpopNRJ <-NA # population en interation avec miel
out$AUCpopDL <-NA #population en interaction avec drone
out$AUCLarveDL <-NA #larve et drone en interaction
out$AUCLarveNRJ <-NA #larve et miel en interaction
out$AUCNRJDL <-NA #larve et drone en interaction


### boucle

i<-1

for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
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
 #auc
    roc1 <- roc(etat_futur[names(fitted(mod101)), "survival_n1"] ~ fitted(mod101)) #MS = Modele Selectionn??
    roc2 <- roc(etat_futur[names(fitted(mod102)), "survival_n1"] ~ fitted(mod102)) #MS = Modele Selectionn??
    roc3 <- roc(etat_futur[names(fitted(mod103)), "survival_n1"] ~ fitted(mod103)) #MS = Modele Selectionn??
    roc4 <- roc(etat_futur[names(fitted(mod104)), "survival_n1"] ~ fitted(mod104)) #MS = Modele Selectionn??
    roc5 <- roc(etat_futur[names(fitted(mod105)), "survival_n1"] ~ fitted(mod105)) #MS = Modele Selectionn??
    roc6 <- roc(etat_futur[names(fitted(mod106)), "survival_n1"] ~ fitted(mod106)) #MS = Modele Selectionn??
    AUCMS1 <- auc(roc1)
    AUCMS2 <- auc(roc2)
    AUCMS3 <- auc(roc3)
    AUCMS4 <- auc(roc4)
    AUCMS5 <- auc(roc5)
    AUCMS6 <- auc(roc6)
    out[i,"AUCpopNRJ"] <- AUCMS2
    out[i,"AUCpopDL"] <- AUCMS3
    out[i,"AUCLarveNRJ"] <- AUCMS4
    out[i,"AUCLarveDL"] <- AUCMS5
    out[i,"AUCNRJDL"] <- AUCMS6

  }
} 

####matrice d'interaction
par(mfrow = c(2,3))

#matrice pop*larve
matAUCpoplarve<- matrix(out$AUCpoplarve,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, matAUCpoplarve, col=rev(cols), main ="AUC interaction pop larve", zlim = c(0.5, 1))
contour(interv,ttime_n, matAUCpoplarve, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")
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

type <- read.table("/Users/solene/Desktop/out graph ect/out_avec_modelestep_et_type.csv", sep=";",dec=".",h=T)
mattype<- matrix(type$type,length(interv),length(ttime_n),byrow=T)
fields::image.plot(interv, ttime_n, mattype, col=rev(cols), main ="Mod??le gard?? 2", nlevel = 5)
contour(interv,ttime_n, mattype, levels = seq(0.5, 1, by = 0.05),add = TRUE, col = "black")


########################################################
############ r??alisation de matrice 3D des interactions
#######################################################

################# 3D au jour 110 (30 jour d'amplitude), au jour 160 (amplitude 20), au jour 200 amplitude 40, au jour 220 (amplitude 150), pour la taille de la population en interaction
#avec les larves, NRJ miel et les larves m??les

###### param??tres d'analyse:

ttime_n=130
# temps ?? pr??dire (intervalle dans le temps)
#interv <- seq(pas, interv_max,pas)
interv=30
###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
i<-1

for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
if(out[i,"prop_survival_n1"]<seuil_survival){
  mod15_130 <- glm(survival_n1 ~ TotalDroneLarvae*TotalLarvae, family = binomial, data = etat_futur)
  mod14_130 <- glm(survival_n1 ~ TotalLarvae*HoneyEnergyStore, family = binomial, data = etat_futur)
}
}

#ttime_n : moment o?? on souhaite regarder pour un indicateur
#ttime_n <- seq(ttime_min,ttime_max,pas)
ttime_n=160
# temps ?? pr??dire (intervalle dans le temps)
#interv <- seq(pas, interv_max,pas)
interv=20
###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
i<-1

for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod14_160 <- glm(survival_n1 ~ HoneyEnergyStore*TotalLarvae, family = binomial, data = etat_futur)
    mod15_160 <- glm(survival_n1 ~ TotalDroneLarvae*TotalLarvae, family = binomial, data = etat_futur)
    
  }
}




#ttime_n : moment o?? on souhaite regarder pour un indicateur
#ttime_n <- seq(ttime_min,ttime_max,pas)
ttime_n=180
# temps ?? pr??dire (intervalle dans le temps)
#interv <- seq(pas, interv_max,pas)
interv=10
###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
i<-1

for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod15_180 <- glm(survival_n1 ~ TotalDroneLarvae*TotalLarvae, family = binomial, data = etat_futur)
    mod14_180 <- glm(survival_n1 ~ TotalLarvae*HoneyEnergyStore, family = binomial, data = etat_futur)
    
  }
}



#ttime_n : moment o?? on souhaite regarder pour un indicateur
#ttime_n <- seq(ttime_min,ttime_max,pas)
ttime_n=220
# temps ?? pr??dire (intervalle dans le temps)
#interv <- seq(pas, interv_max,pas)
interv=140
###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
i<-1

for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod14_220 <- glm(survival_n1 ~ HoneyEnergyStore*TotalLarvae, family = binomial, data = etat_futur)
    mod15_220 <- glm(survival_n1 ~ TotalDroneLarvae*TotalLarvae, family = binomial, data = etat_futur)
    
  }
}

#######r??alisation des courbes 3D

par(mfrow = c(1,1))
##3D
#POP*LARVE
visreg2d(mod11_110, x="TotalPopSize", y="TotalLarvae", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="Nombre de Larves", main = "3D Population * Larves ?? 120 jours pr??vision ?? 30")

visreg2d(mod11_160, x="TotalPopSize", y="TotalLarvae", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="Nombre de Larves", main = "3D Population * Larves ?? 160 jours pr??vision ?? 20")

visreg2d(mod11_200, x="TotalPopSize", y="TotalLarvae", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="Nombre de Larves", main = "3D Population * Larves ?? 200 jours pr??vision ?? 40")

visreg2d(mod11_220, x="TotalPopSize", y="TotalLarvae", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="Nombre de Larves", main = "3D Population * Larves ?? 220 jours pr??vision ?? 150")

par(mfrow = c(2,2))
#POP*NRJ
visreg2d(mod12_110, x="TotalPopSize", y="HoneyEnergyStore", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="NRJ (miel)", main = "3D Population * HoneyEnergyStore ?? 120 jours pr??vision ?? 30")

visreg2d(mod12_160, x="TotalPopSize", y="HoneyEnergyStore", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="NRJ (miel)", main = "3D Population * HoneyEnergyStore ?? 160 jours pr??vision ?? 20")

visreg2d(mod12_200, x="TotalPopSize", y="HoneyEnergyStore", plot.type='persp',scale='response', col="white", border="black", 
         theta=180, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="NRJ (miel)", main = "3D Population * HoneyEnergyStore ?? 200 jours pr??vision ?? 40")

visreg2d(mod12_220, x="TotalPopSize", y="HoneyEnergyStore", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="NRJ (miel)", main = "3D Population * HoneyEnergyStore ?? 220 jours pr??vision ?? 150")

par(mfrow = c(2,2))
#POP*DL
visreg2d(mod13_110, x="TotalPopSize", y="TotalDroneLarvae", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="Nombre de Larves m??le", main = "3D Population * TotalDroneLarvae ?? 120 jours pr??vision ?? 30")

visreg2d(mod13_160, x="TotalPopSize", y="TotalDroneLarvae", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="Nombre de Larves m??le", main = "3D Population * TotalDroneLarvae ?? 160 jours pr??vision ?? 20")

visreg2d(mod13_200, x="TotalPopSize", y="TotalDroneLarvae", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="Nombre de Larves m??le", main = "3D Population * TotalDroneLarvae ?? 200 jours pr??vision ?? 40")

visreg2d(mod13_220, x="TotalPopSize", y="TotalDroneLarvae", plot.type='persp',scale='response', col="white", border="black", 
         theta=60, phi=10, r=10, zlab="Survie", xlab="Taille de la population", ylab="Nombre de Larves m??le", main = "3D Population * TotalDroneLarvae ?? 220 jours pr??vision ?? 150")



#################################################################
#######################courbe r??capitulative jusqu'?? 60 jours de pr??dictions
################### 
#le fichier avec les mod??les d??j?? fait :

out <- read.table("/Users/solene/Desktop/out graph ect/fichier excel out/Out_graph_courbe_90.csv", sep=";",dec=".",h=T)

### ?? 90 jours
# ttime_min : jour o?? commencent les analyses
ttime_min<-100 
#ttime_max : jour o?? finissent les analyses
ttime_max<-130
# combien de temps plus tard on regarde 
pas <- 10

seuil_survival <- 0.99 #Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse

####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 60

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- 100
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
head(out)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i<-seq(1,i_max,1)
out[i,"nb_survival"] <- sum(etat_futur$survival)
out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
out$AUCpoidstot_90 <- NA
out$AUClarve_90 <- NA
out$AUCLD_90 <- NA
out$AUCpop_90 <- NA
out$AUCmiel_90 <- NA
out$AUClarveNRJ_90 <- NA
out$AUClarveLD_90 <- NA
out$AUCAFF_90 <- NA
out$AUCLSP_90 <- NA
out$AUCTotevend_90 <- NA

### boucle

#i<-4
i <-1
for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod1_90 <- glm(survival_n1 ~ TotalWeight, family = binomial, data = etat_futur)
    mod2_90 <- glm(survival_n1 ~ TotalLarvae, family = binomial, data = etat_futur)
    #mod3_90 <- glm(survival_n1 ~ TotalDroneLarvae, family = binomial, data = etat_futur)
    mod4_90 <- glm(survival_n1 ~ TotalPopSize, family = binomial, data = etat_futur)
    mod5_90 <- glm(survival_n1 ~ HoneyEnergyStore, family = binomial, data = etat_futur)
    mod6_90 <- glm(survival_n1 ~ (TotalLarvae+HoneyEnergyStore)+(TotalLarvae:HoneyEnergyStore), family = binomial, data = etat_futur)
    #mod7_90 <- glm(survival_n1 ~ TotalLarvae*TotalDroneLarvae, family = binomial, data = etat_futur)
    mod8_90 <- glm(survival_n1 ~ AFF, family = binomial, data = etat_futur)
    mod9_90 <- glm(survival_n1 ~ LSP, family = binomial, data = etat_futur)
    mod10_90 <- glm(survival_n1 ~ TotalEventsToday, family = binomial, data = etat_futur)
    
    roc1_90 <- roc(etat_futur[names(fitted(mod1_90)), "survival_n1"] ~ fitted(mod1_90))
    roc2_90 <- roc(etat_futur[names(fitted(mod2_90)), "survival_n1"] ~ fitted(mod2_90))
    roc3_90 <- roc(etat_futur[names(fitted(mod3_90)), "survival_n1"] ~ fitted(mod3_90))
    roc4_90 <- roc(etat_futur[names(fitted(mod4_90)), "survival_n1"] ~ fitted(mod4_90))
    roc5_90 <- roc(etat_futur[names(fitted(mod5_90)), "survival_n1"] ~ fitted(mod5_90))
    roc6_90 <- roc(etat_futur[names(fitted(mod6_90)), "survival_n1"] ~ fitted(mod6_90))
    roc7_90 <- roc(etat_futur[names(fitted(mod7_90)), "survival_n1"] ~ fitted(mod7_90))
    roc8_90 <- roc(etat_futur[names(fitted(mod8_90)), "survival_n1"] ~ fitted(mod8_90))
    roc9_90 <- roc(etat_futur[names(fitted(mod9_90)), "survival_n1"] ~ fitted(mod9_90))
    roc10_90 <- roc(etat_futur[names(fitted(mod10_90)), "survival_n1"] ~ fitted(mod10_90))
    
    AUC1_90 <- auc(roc1_90)
    out[i,"AUCpoidstot_90"] <- AUC1_90
    AUC2_90 <- auc(roc2_90)
    out[i,"AUClarve_90"] <- AUC2_90
    AUC3_90 <- auc(roc3_90)
    out[i,"AUCLD_90"] <- AUC3_90
    AUC4_90 <- auc(roc4_90)
    out[i,"AUCpop_90"] <- AUC4_90
    AUC5_90 <- auc(roc5_90)
    out[i,"AUCmiel_90"] <- AUC5_90   
    AUC6_90 <- auc(roc6_90)
    out[i,"AUClarveNRJ_90"] <- AUC6_90
    AUC7_90 <- auc(roc7_90)
    out[i,"AUClarveLD_90"] <- AUC7_90
    AUC8_90 <- auc(roc8_90)
    out[i,"AUCAFF_90"] <- AUC8_90
    AUC9_90 <- auc(roc9_90)
    out[i,"AUCLSP_90"] <- AUC9_90
    AUC10_90 <- auc(roc10_90)
    out[i,"AUCTotevend_90"] <- AUC10_90
  }
}

par(mfrow = c(1,1))


plot(out$interv,out$AUCpoidstot_90,lwd=3, lty=3, type='l',las=1,col="red", ylab="Pouvoir pr??dictif (AUC)", xlab="Port??e pr??dictive (en jours)", main="a) D??but de la p??riode du colza", ylim=c(0.5,1))
lines(out$interv,out$AUClarve_90, col="dark green", lwd=3)
#lines(out$interv,out$AUCLD_90, col="dark blue", lwd=3)
lines(out$interv,out$AUCpop_90, col="blue", lwd=3)
lines(out$interv,out$AUCmiel_90, col="green", lwd=3)
lines(out$interv,out$AUClarveNRJ_90, col="purple", lwd=3, lty=5)
#lines(out$interv,out$AUClarveLD_90, col="black", lwd=3, lty=5)
lines(out$interv,out$AUCAFF_90, col="yellow", lwd=3, lty=4)
lines(out$interv,out$AUCLSP_90, col="orange", lwd=3, lty=4)
lines(out$interv,out$AUCTotevend_90, col="pink", lwd=3, lty=4)

legend("topright",inset =c(0.12,0),xjust=0,yjust=0, legend="Poids de la colonie",col="red", lty=3,lwd=3, cex=0.7, title ="Poids de la colonie", box.lty=0)
legend("topright", inset =c(0.04,0.08), box.lty=0,xjust=0,yjust=0, legend=c("Couvain femelle", "Couvain m??le", "Taille population", "Energie miel"),col=c(" dark green", "dark blue","blue", "green"), lty=1,lwd=3, cex=0.7, title="Traits d??mographiques s??par??ment")
legend("topright", box.lty=0,inset =c(0.01,0.25),xjust=0,yjust=0,legend=c("Couvain*Energie miel", "Couvain femelle*Couvain m??le"), col=c("purple","black"), lty=5,lwd=3, cex=0.7, title="Traits d??mographiques en interaction")
legend("topright",box.lty=0,inset =c(0.07,0.35),xjust=0,yjust=0,legend=c("AFF", "LSP", "Activit?? de vols"), col=c("yellow","orange", "pink") ,lty=4,lwd=3, cex=0.7, title ="Traits de vie individuels")

########en anglais
plot(out$interv,out$AUCpoidstot_90,lwd=3, lty=3, type='l', col="red", ylab="Predictive power (AUC)", xlab="Prediction timespan (in days)",main="Predictive power and timespan for different colony 
     collapse candidate predictors (starting at day 180 of the year)" , ylim=c(0.5,1))
lines(out$interv,out$AUClarve_90, col="dark green", lwd=3)
lines(out$interv,out$AUCLD_90, col="dark blue", lwd=3)
lines(out$interv,out$AUCpop_90, col="blue", lwd=3)
lines(out$interv,out$AUCmiel_90, col="green", lwd=3)
lines(out$interv,out$AUClarveNRJ_90, col="purple", lwd=3, lty=5)
lines(out$interv,out$AUClarveLD_90, col="black", lwd=3, lty=5)
lines(out$interv,out$AUCAFF_90, col="yellow", lwd=3, lty=4)
lines(out$interv,out$AUCLSP_90, col="orange", lwd=3, lty=4)
lines(out$interv,out$AUCTotevend_90, col="pink", lwd=3, lty=4)
legend("topright",inset =c(0.12,0),xjust=0,yjust=0, legend="Colony mass",col="red", lty=3,lwd=3, cex=0.7, title ="Colony mass", box.lty=0)
legend("topright", inset =c(0.04,0.08), box.lty=0,xjust=0,yjust=0, legend=c("Worker Brood", "Drone Brood", "Adult population", "Honey reserves"),col=c(" dark green", "dark blue","blue", "green"), lty=1,lwd=3, cex=0.7, title="Seperate Colony level traits")
legend("topright", box.lty=0,inset =c(0.01,0.25),xjust=0,yjust=0,legend=c("Worker Brood*Honey reserve", "Worker Brood*Drone Brood"), col=c("purple","black"), lty=5,lwd=3, cex=0.7, title="Interacting Colony traits")
legend("topright",box.lty=0,inset =c(0.07,0.35),xjust=0,yjust=0,legend=c("AFF", "LSP", "Daily flight activity"), col=c("yellow","orange", "pink") ,lty=4,lwd=3, cex=0.7, title ="Individual level traits")


 
### ?? 120 jours
# ttime_min : jour o?? commencent les analyses
ttime_min<-100 
#ttime_max : jour o?? finissent les analyses
ttime_max<-130
# combien de temps plus tard on regarde 
pas <- 10

seuil_survival <- 0.99 #Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse

####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 150

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- 120
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
head(out)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i<-seq(1,i_max,1)
out[i,"nb_survival"] <- sum(etat_futur$survival)
out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
out$AUCpoidstot_120 <- NA
out$AUClarve_120 <- NA
out$AUCLD_120 <- NA
out$AUCpop_120 <- NA
out$AUCmiel_120 <- NA
out$AUClarveNRJ_120 <- NA
out$AUClarveLD_120 <- NA
out$AUCAFF_120 <- NA
out$AUCLSP_120 <- NA
out$AUCTotevend_120 <- NA

### boucle

#i<-4
i <-1
for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod1_120 <- glm(survival_n1 ~ TotalWeight, family = binomial, data = etat_futur)
    mod2_120 <- glm(survival_n1 ~ TotalLarvae, family = binomial, data = etat_futur)
    mod3_120 <- glm(survival_n1 ~ TotalDroneLarvae, family = binomial, data = etat_futur)
    mod4_120 <- glm(survival_n1 ~ TotalPopSize, family = binomial, data = etat_futur)
    mod5_120 <- glm(survival_n1 ~ HoneyEnergyStore, family = binomial, data = etat_futur)
    mod6_120 <- glm(survival_n1 ~ TotalLarvae*HoneyEnergyStore, family = binomial, data = etat_futur)
    mod7_120 <- glm(survival_n1 ~ TotalLarvae*TotalDroneLarvae, family = binomial, data = etat_futur)
    mod8_120 <- glm(survival_n1 ~ AFF, family = binomial, data = etat_futur)
    mod9_120 <- glm(survival_n1 ~ LSP, family = binomial, data = etat_futur)
    mod10_120 <- glm(survival_n1 ~ TotalEventsToday, family = binomial, data = etat_futur)
    
    roc1_120 <- roc(etat_futur[names(fitted(mod1_120)), "survival_n1"] ~ fitted(mod1_120))
    roc2_120 <- roc(etat_futur[names(fitted(mod2_120)), "survival_n1"] ~ fitted(mod2_120))
    roc3_120 <- roc(etat_futur[names(fitted(mod3_120)), "survival_n1"] ~ fitted(mod3_120))
    roc4_120 <- roc(etat_futur[names(fitted(mod4_120)), "survival_n1"] ~ fitted(mod4_120))
    roc5_120 <- roc(etat_futur[names(fitted(mod5_120)), "survival_n1"] ~ fitted(mod5_120))
    roc6_120 <- roc(etat_futur[names(fitted(mod6_120)), "survival_n1"] ~ fitted(mod6_120))
    roc7_120 <- roc(etat_futur[names(fitted(mod7_120)), "survival_n1"] ~ fitted(mod7_120))
    roc8_120 <- roc(etat_futur[names(fitted(mod8_120)), "survival_n1"] ~ fitted(mod8_120))
    roc9_120 <- roc(etat_futur[names(fitted(mod9_120)), "survival_n1"] ~ fitted(mod9_120))
    roc10_120 <- roc(etat_futur[names(fitted(mod10_120)), "survival_n1"] ~ fitted(mod10_120))
    
    AUC1_120 <- auc(roc1_120)
    out[i,"AUCpoidstot_120"] <- AUC1_120
    AUC2_120 <- auc(roc2_120)
    out[i,"AUClarve_120"] <- AUC2_120
    AUC3_120 <- auc(roc3_120)
    out[i,"AUCLD_120"] <- AUC3_120
    AUC4_120 <- auc(roc4_120)
    out[i,"AUCpop_120"] <- AUC4_120
    AUC5_120 <- auc(roc5_120)
    out[i,"AUCmiel_120"] <- AUC5_120   
    AUC6_120 <- auc(roc6_120)
    out[i,"AUClarveNRJ_120"] <- AUC6_120
    AUC7_120 <- auc(roc7_120)
    out[i,"AUClarveLD_120"] <- AUC7_120
    AUC8_120 <- auc(roc8_120)
    out[i,"AUCAFF_120"] <- AUC8_120
    AUC9_120 <- auc(roc9_120)
    out[i,"AUCLSP_120"] <- AUC9_120
    AUC10_120 <- auc(roc10_120)
    out[i,"AUCTotevend_120"] <- AUC10_120
  }
}

plot(out$interv,out$AUCpoidstot_120,lwd=3, lty=3, type='l', col="red", ylab="Puissance (AUC)", xlab="Pr??diction (en jours)",main="Puissance et port??e 
des diff??rents indicateurs au jour d'observation 120 (avant Disette)", ylim=c(0.5,1))
lines(out$interv,out$AUClarve_120, col="dark green", lwd=3)
lines(out$interv,out$AUCLD_120, col="dark blue", lwd=3)
lines(out$interv,out$AUCpop_120, col="blue", lwd=3)
lines(out$interv,out$AUCmiel_120, col="green", lwd=3)
lines(out$interv,out$AUClarveNRJ_120, col="purple", lwd=3, lty=5)
lines(out$interv,out$AUClarveLD_120, col="black", lwd=3, lty=5)
lines(out$interv,out$AUCAFF_120, col="yellow", lwd=3, lty=4)
lines(out$interv,out$AUCLSP_120, col="orange", lwd=3, lty=4)
lines(out$interv,out$AUCTotevend_120, col="pink", lwd=3, lty=4)

legend("topright",inset =c(0.12,0),xjust=0,yjust=0, legend="Poids de la colonie",col="red", lty=3,lwd=3, cex=0.7, title ="Poids de la colonie", box.lty=0)
legend("topright", inset =c(0.04,0.08), box.lty=0,xjust=0,yjust=0, legend=c("Couvain femelle", "Couvain m??le", "Taille population", "Energie miel"),col=c(" dark green", "dark blue","blue", "green"), lty=1,lwd=3, cex=0.7, title="Traits d??mographiques s??par??ment")
legend("topright", box.lty=0,inset =c(0.01,0.25),xjust=0,yjust=0,legend=c("Couvain*Quantit?? miel stock??ee", "Couvain femelle*Couvain m??le"), col=c("purple","black"), lty=5,lwd=3, cex=0.7, title="Traits d??mographiques en interaction")
legend("topright",box.lty=0,inset =c(0.07,0.35),xjust=0,yjust=0,legend=c("AFF", "LSP", "Activit?? de vols"), col=c("yellow","orange", "pink") ,lty=4,lwd=3, cex=0.7, title ="Traits de vie individuels")


### ?? 180 jours
# ttime_min : jour o?? commencent les analyses
ttime_min<-100 
#ttime_max : jour o?? finissent les analyses
ttime_max<-130
# combien de temps plus tard on regarde 
pas <- 10

seuil_survival <- 0.99 #Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse

####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 60

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- 180
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
head(out)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i<-seq(1,i_max,1)
out[i,"nb_survival"] <- sum(etat_futur$survival)
out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
out$AUCpoidstot_180 <- NA
out$AUClarve_180 <- NA
out$AUCLD_180 <- NA
out$AUCpop_180 <- NA
out$AUCmiel_180 <- NA
out$AUClarveNRJ_180 <- NA
out$AUClarveLD_180 <- NA
out$AUCAFF_180 <- NA
out$AUCLSP_180 <- NA
out$AUCTotevend_180 <- NA

### boucle

#i<-4
i <-1
for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod1_180 <- glm(survival_n1 ~ TotalWeight, family = binomial, data = etat_futur)
    mod2_180 <- glm(survival_n1 ~ TotalLarvae, family = binomial, data = etat_futur)
    mod3_180 <- glm(survival_n1 ~ TotalDroneLarvae, family = binomial, data = etat_futur)
    mod4_180 <- glm(survival_n1 ~ TotalPopSize, family = binomial, data = etat_futur)
    mod5_180 <- glm(survival_n1 ~ HoneyEnergyStore, family = binomial, data = etat_futur)
    mod6_180 <- glm(survival_n1 ~ TotalLarvae*HoneyEnergyStore, family = binomial, data = etat_futur)
    mod7_180 <- glm(survival_n1 ~ TotalLarvae*TotalDroneLarvae, family = binomial, data = etat_futur)
    mod8_180 <- glm(survival_n1 ~ AFF, family = binomial, data = etat_futur)
    mod9_180 <- glm(survival_n1 ~ LSP, family = binomial, data = etat_futur)
    mod10_180 <- glm(survival_n1 ~ TotalEventsToday, family = binomial, data = etat_futur)
    
    roc1_180 <- roc(etat_futur[names(fitted(mod1_180)), "survival_n1"] ~ fitted(mod1_180))
    roc2_180 <- roc(etat_futur[names(fitted(mod2_180)), "survival_n1"] ~ fitted(mod2_180))
    roc3_180 <- roc(etat_futur[names(fitted(mod3_180)), "survival_n1"] ~ fitted(mod3_180))
    roc4_180 <- roc(etat_futur[names(fitted(mod4_180)), "survival_n1"] ~ fitted(mod4_180))
    roc5_180 <- roc(etat_futur[names(fitted(mod5_180)), "survival_n1"] ~ fitted(mod5_180))
    roc6_180 <- roc(etat_futur[names(fitted(mod6_180)), "survival_n1"] ~ fitted(mod6_180))
    roc7_180 <- roc(etat_futur[names(fitted(mod7_180)), "survival_n1"] ~ fitted(mod7_180))
    roc8_180 <- roc(etat_futur[names(fitted(mod8_180)), "survival_n1"] ~ fitted(mod8_180))
    roc9_180 <- roc(etat_futur[names(fitted(mod9_180)), "survival_n1"] ~ fitted(mod9_180))
    roc10_180 <- roc(etat_futur[names(fitted(mod10_180)), "survival_n1"] ~ fitted(mod10_180))
    
    AUC1_180 <- auc(roc1_180)
    out[i,"AUCpoidstot_180"] <- AUC1_180
    AUC2_180 <- auc(roc2_180)
    out[i,"AUClarve_180"] <- AUC2_180
    AUC3_180 <- auc(roc3_180)
    out[i,"AUCLD_180"] <- AUC3_180
    AUC4_180 <- auc(roc4_180)
    out[i,"AUCpop_180"] <- AUC4_180
    AUC5_180 <- auc(roc5_180)
    out[i,"AUCmiel_180"] <- AUC5_180   
    AUC6_180 <- auc(roc6_180)
    out[i,"AUClarveNRJ_180"] <- AUC6_180
    AUC7_180 <- auc(roc7_180)
    out[i,"AUClarveLD_180"] <- AUC7_180
    AUC8_180 <- auc(roc8_180)
    out[i,"AUCAFF_180"] <- AUC8_180
    AUC9_180 <- auc(roc9_180)
    out[i,"AUCLSP_180"] <- AUC9_180
    AUC10_180 <- auc(roc10_180)
    out[i,"AUCTotevend_180"] <- AUC10_180
  }
}

plot(out$interv,out$AUCpoidstot_180,lwd=3, lty=3, type='l', col="red", ylab="Pouvoir pr??dictif (AUC)", las=1, xlab="Port??e pr??dictive (en jours)", main="b) D??but de la p??riode tournesol/ma??s", ylim=c(0.5,1))
lines(out$interv,out$AUClarve_180, col="dark green", lwd=3)
lines(out$interv,out$AUCLD_180, col="dark blue", lwd=3)
lines(out$interv,out$AUCpop_180, col="blue", lwd=3)
lines(out$interv,out$AUCmiel_180, col="green", lwd=3)
lines(out$interv,out$AUClarveNRJ_180, col="purple", lwd=3, lty=5)
lines(out$interv,out$AUClarveLD_180, col="black", lwd=3, lty=5)
lines(out$interv,out$AUCAFF_180, col="yellow", lwd=3, lty=4)
lines(out$interv,out$AUCLSP_180, col="orange", lwd=3, lty=4)
lines(out$interv,out$AUCTotevend_180, col="pink", lwd=3, lty=4)

legend("topright",inset =c(0.12,0),xjust=0,yjust=0, legend="Poids de la colonie",col="red", lty=3,lwd=3, cex=0.7, title ="Poids de la colonie", box.lty=0)
legend("topright", inset =c(0.04,0.08), box.lty=0,xjust=0,yjust=0, legend=c("Couvain femelle", "Couvain m??le", "Taille population", "Energie miel"),col=c(" dark green", "dark blue","blue", "green"), lty=1,lwd=3, cex=0.7, title="Traits d??mographiques s??par??ment")
legend("topright", box.lty=0,inset =c(0.01,0.25),xjust=0,yjust=0,legend=c("Couvain*Energie miel", "Couvain femelle*Couvain m??le"), col=c("purple","black"), lty=5,lwd=3, cex=0.7, title="Traits d??mographiques en interaction")
legend("topright",box.lty=0,inset =c(0.07,0.35),xjust=0,yjust=0,legend=c("AFF", "LSP", "Activit?? de vols"), col=c("yellow","orange", "pink") ,lty=4,lwd=3, cex=0.7, title ="Traits de vie individuels")



### ?? 210 jours
# ttime_min : jour o?? commencent les analyses
ttime_min<-100 
#ttime_max : jour o?? finissent les analyses
ttime_max<-130
# combien de temps plus tard on regarde 
pas <- 10

seuil_survival <- 0.99 #Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse

####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 150

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- 210
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
head(out)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i<-seq(1,i_max,1)
out[i,"nb_survival"] <- sum(etat_futur$survival)
out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
out$AUCpoidstot_210 <- NA
out$AUClarve_210 <- NA
out$AUCLD_210 <- NA
out$AUCpop_210 <- NA
out$AUCmiel_210 <- NA
out$AUClarveNRJ_210 <- NA
out$AUClarveLD_210 <- NA
out$AUCAFF_210 <- NA
out$AUCLSP_210 <- NA
out$AUCTotevend_210 <- NA

### boucle

#i<-4
i <-1
for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod1_210 <- glm(survival_n1 ~ TotalWeight, family = binomial, data = etat_futur)
    mod2_210 <- glm(survival_n1 ~ TotalLarvae, family = binomial, data = etat_futur)
    mod3_210 <- glm(survival_n1 ~ TotalDroneLarvae, family = binomial, data = etat_futur)
    mod4_210 <- glm(survival_n1 ~ TotalPopSize, family = binomial, data = etat_futur)
    mod5_210 <- glm(survival_n1 ~ HoneyEnergyStore, family = binomial, data = etat_futur)
    mod6_210 <- glm(survival_n1 ~ TotalLarvae*HoneyEnergyStore, family = binomial, data = etat_futur)
    mod7_210 <- glm(survival_n1 ~ TotalLarvae*TotalDroneLarvae, family = binomial, data = etat_futur)
    mod8_210 <- glm(survival_n1 ~ AFF, family = binomial, data = etat_futur)
    mod9_210 <- glm(survival_n1 ~ LSP, family = binomial, data = etat_futur)
    mod10_210 <- glm(survival_n1 ~ TotalEventsToday, family = binomial, data = etat_futur)
    
    roc1_210 <- roc(etat_futur[names(fitted(mod1_210)), "survival_n1"] ~ fitted(mod1_210))
    roc2_210 <- roc(etat_futur[names(fitted(mod2_210)), "survival_n1"] ~ fitted(mod2_210))
    roc3_210 <- roc(etat_futur[names(fitted(mod3_210)), "survival_n1"] ~ fitted(mod3_210))
    roc4_210 <- roc(etat_futur[names(fitted(mod4_210)), "survival_n1"] ~ fitted(mod4_210))
    roc5_210 <- roc(etat_futur[names(fitted(mod5_210)), "survival_n1"] ~ fitted(mod5_210))
    roc6_210 <- roc(etat_futur[names(fitted(mod6_210)), "survival_n1"] ~ fitted(mod6_210))
    roc7_210 <- roc(etat_futur[names(fitted(mod7_210)), "survival_n1"] ~ fitted(mod7_210))
    roc8_210 <- roc(etat_futur[names(fitted(mod8_210)), "survival_n1"] ~ fitted(mod8_210))
    roc9_210 <- roc(etat_futur[names(fitted(mod9_210)), "survival_n1"] ~ fitted(mod9_210))
    roc10_210 <- roc(etat_futur[names(fitted(mod10_210)), "survival_n1"] ~ fitted(mod10_210))
    
    AUC1_210 <- auc(roc1_210)
    out[i,"AUCpoidstot_210"] <- AUC1_210
    AUC2_210 <- auc(roc2_210)
    out[i,"AUClarve_210"] <- AUC2_210
    AUC3_210 <- auc(roc3_210)
    out[i,"AUCLD_210"] <- AUC3_210
    AUC4_210 <- auc(roc4_210)
    out[i,"AUCpop_210"] <- AUC4_210
    AUC5_210 <- auc(roc5_210)
    out[i,"AUCmiel_210"] <- AUC5_210   
    AUC6_210 <- auc(roc6_210)
    out[i,"AUClarveNRJ_210"] <- AUC6_210
    AUC7_210 <- auc(roc7_210)
    out[i,"AUClarveLD_210"] <- AUC7_210
    AUC8_210 <- auc(roc8_210)
    out[i,"AUCAFF_210"] <- AUC8_210
    AUC9_210 <- auc(roc9_210)
    out[i,"AUCLSP_210"] <- AUC9_210
    AUC10_210 <- auc(roc10_210)
    out[i,"AUCTotevend_210"] <- AUC10_210
  }
}

plot(out$interv,out$AUCpoidstot_210,lwd=3, lty=3, type='l', col="red", ylab="Puissance (AUC)", xlab="Pr??diction (en jours)",main="Puissance et port??e 
     des diff??rents indicateurs au jour d'observation 210 (avant Ma??s)", ylim=c(0.5,1))
lines(out$interv,out$AUClarve_210, col="dark green", lwd=3)
lines(out$interv,out$AUCLD_210, col="dark blue", lwd=3)
lines(out$interv,out$AUCpop_210, col="blue", lwd=3)
lines(out$interv,out$AUCmiel_210, col="green", lwd=3)
lines(out$interv,out$AUClarveNRJ_210, col="purple", lwd=3, lty=5)
lines(out$interv,out$AUClarveLD_210, col="black", lwd=3, lty=5)
lines(out$interv,out$AUCAFF_210, col="yellow", lwd=3, lty=4)
lines(out$interv,out$AUCLSP_210, col="orange", lwd=3, lty=4)
lines(out$interv,out$AUCTotevend_210, col="pink", lwd=3, lty=4)

legend("topright",inset =c(0.12,0),xjust=0,yjust=0, legend="Poids de la colonie",col="red", lty=3,lwd=3, cex=0.7, title ="Poids de la colonie", box.lty=0)
legend("topright", inset =c(0.04,0.08), box.lty=0,xjust=0,yjust=0, legend=c("Couvain ouvri??re", "Couvain m??le", "Taille population", "Energie miel"),col=c(" dark green", "dark blue","blue", "green"), lty=1,lwd=3, cex=0.7, title="Traits d??mographiques s??par??ment")
legend("topright", box.lty=0,inset =c(0.01,0.25),xjust=0,yjust=0,legend=c("Couvain ouvri??re*Energie miel", "Couvain ouvri??re*Couvain m??le"), col=c("purple","black"), lty=5,lwd=3, cex=0.7, title="Traits d??mographiques en interaction")
legend("topright",box.lty=0,inset =c(0.07,0.35),xjust=0,yjust=0,legend=c("AFF", "LSP", "Activit?? de vols"), col=c("yellow","orange", "pink") ,lty=4,lwd=3, cex=0.7, title ="Traits de vie individuels")



### ?? 390 jours
# ttime_min : jour o?? commencent les analyses
ttime_min<-100 
#ttime_max : jour o?? finissent les analyses
ttime_max<-130
# combien de temps plus tard on regarde 
pas <- 10

seuil_survival <- 0.99 #Seuil max de survie entre "ttime_n" et "ttime_n+interv" pour faire l'analyse

####### param??tre calcul??
#interv_max : le plus grand intervalle possible
interv_max <- 150

#ttime_n : moment o?? on souhaite regarder pour un indicateur
ttime_n <- 390
# temps ?? pr??dire (intervalle dans le temps)
interv <- seq(pas, interv_max,pas)

###### Initialisation de la boucle: pr??paration des objets de sortie

#out : creation de la grille de sortie avec ttime_n et interv
out <- expand.grid(ttime_n=ttime_n,interv=interv)
head(out)
#i_max : nombre de lignes dans out = nb total d'analyses = nb de cellumes dans la matrice finale
i_max <- nrow(out) 
out$i<-seq(1,i_max,1)
out[i,"nb_survival"] <- sum(etat_futur$survival)
out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
out$AUCpoidstot_390 <- NA
out$AUClarve_390 <- NA
out$AUCLD_390 <- NA
out$AUCpop_390 <- NA
out$AUCmiel_390 <- NA
out$AUClarveNRJ_390 <- NA
out$AUClarveLD_390 <- NA
out$AUCAFF_390 <- NA
out$AUCLSP_390 <- NA
out$AUCTotevend_390 <- NA

### boucle

#i<-4
i <-1
for(i in 1:i_max) {
  time_i <- out[i,"ttime_n"]
  interv_i <- out[i,"interv"]
  
  etat_i <- subset(simulRes,ttime==time_i&survival==1) #donn??e pr??sente en selectionnant que les colonies qui survivent
  futur_i <- subset(simulRes,ttime==time_i+interv_i)[,c("param","survival")] # donn??e future 
  colnames(futur_i)[2]<-"survival_n1" # changer le nom de la colonne survival dans futur_i
  
  etat_futur<-merge(etat_i,futur_i,by="param") #importation des donn??es survie future dans un m??me tableau
  
  out[i,"nb_survival"] <- sum(etat_futur$survival)
  out[i,"nb_survival_n1"] <- sum(etat_futur$survival_n1)
  out[i,"prop_survival_n1"] <- out[i,"nb_survival_n1"]/out[i,"nb_survival"]
  
  if(out[i,"prop_survival_n1"]<seuil_survival){
    mod1_390 <- glm(survival_n1 ~ TotalWeight, family = binomial, data = etat_futur)
    mod2_390 <- glm(survival_n1 ~ TotalLarvae, family = binomial, data = etat_futur)
    mod3_390 <- glm(survival_n1 ~ TotalDroneLarvae, family = binomial, data = etat_futur)
    mod4_390 <- glm(survival_n1 ~ TotalPopSize, family = binomial, data = etat_futur)
    mod5_390 <- glm(survival_n1 ~ HoneyEnergyStore, family = binomial, data = etat_futur)
    mod6_390 <- glm(survival_n1 ~ TotalLarvae*HoneyEnergyStore, family = binomial, data = etat_futur)
    mod7_390 <- glm(survival_n1 ~ TotalLarvae*TotalDroneLarvae, family = binomial, data = etat_futur)
    mod8_390 <- glm(survival_n1 ~ AFF, family = binomial, data = etat_futur)
    mod9_390 <- glm(survival_n1 ~ LSP, family = binomial, data = etat_futur)
    mod10_390 <- glm(survival_n1 ~ TotalEventsToday, family = binomial, data = etat_futur)
    
    roc1_390 <- roc(etat_futur[names(fitted(mod1_390)), "survival_n1"] ~ fitted(mod1_390))
    roc2_390 <- roc(etat_futur[names(fitted(mod2_390)), "survival_n1"] ~ fitted(mod2_390))
    roc3_390 <- roc(etat_futur[names(fitted(mod3_390)), "survival_n1"] ~ fitted(mod3_390))
    roc4_390 <- roc(etat_futur[names(fitted(mod4_390)), "survival_n1"] ~ fitted(mod4_390))
    roc5_390 <- roc(etat_futur[names(fitted(mod5_390)), "survival_n1"] ~ fitted(mod5_390))
    roc6_390 <- roc(etat_futur[names(fitted(mod6_390)), "survival_n1"] ~ fitted(mod6_390))
    roc7_390 <- roc(etat_futur[names(fitted(mod7_390)), "survival_n1"] ~ fitted(mod7_390))
    roc8_390 <- roc(etat_futur[names(fitted(mod8_390)), "survival_n1"] ~ fitted(mod8_390))
    roc9_390 <- roc(etat_futur[names(fitted(mod9_390)), "survival_n1"] ~ fitted(mod9_390))
    roc10_390 <- roc(etat_futur[names(fitted(mod10_390)), "survival_n1"] ~ fitted(mod10_390))
    
    AUC1_390 <- auc(roc1_390)
    out[i,"AUCpoidstot_390"] <- AUC1_390
    AUC2_390 <- auc(roc2_390)
    out[i,"AUClarve_390"] <- AUC2_390
    AUC3_390 <- auc(roc3_390)
    out[i,"AUCLD_390"] <- AUC3_390
    AUC4_390 <- auc(roc4_390)
    out[i,"AUCpop_390"] <- AUC4_390
    AUC5_390 <- auc(roc5_390)
    out[i,"AUCmiel_390"] <- AUC5_390   
    AUC6_390 <- auc(roc6_390)
    out[i,"AUClarveNRJ_390"] <- AUC6_390
    AUC7_390 <- auc(roc7_390)
    out[i,"AUClarveLD_390"] <- AUC7_390
    AUC8_390 <- auc(roc8_390)
    out[i,"AUCAFF_390"] <- AUC8_390
    AUC9_390 <- auc(roc9_390)
    out[i,"AUCLSP_390"] <- AUC9_390
    AUC10_390 <- auc(roc10_390)
    out[i,"AUCTotevend_390"] <- AUC10_390
  }
}

plot(out$interv,out$AUCpoidstot_390,lwd=3, lty=3, type='l', col="red", ylab="Puissance (AUC)", xlab="Pr??diction (en jours)",main="Puissance et port??e 
     des diff??rents indicateurs au jour d'observation ?? 390 (sortie Hivernage)", ylim=c(0.5,1))
lines(out$interv,out$AUClarve_390, col="dark green", lwd=3)
lines(out$interv,out$AUCLD_390, col="dark blue", lwd=3)
lines(out$interv,out$AUCpop_390, col="blue", lwd=3)
lines(out$interv,out$AUCmiel_390, col="green", lwd=3)
lines(out$interv,out$AUClarveNRJ_390, col="purple", lwd=3, lty=5)
lines(out$interv,out$AUClarveLD_390, col="black", lwd=3, lty=5)
lines(out$interv,out$AUCAFF_390, col="yellow", lwd=3, lty=4)
lines(out$interv,out$AUCLSP_390, col="orange", lwd=3, lty=4)
lines(out$interv,out$AUCTotevend_390, col="pink", lwd=3, lty=4)

legend("topright",inset =c(0.12,0),xjust=0,yjust=0, legend="Poids de la colonie",col="red", lty=3,lwd=3, cex=0.7, title ="Poids de la colonie", box.lty=0)
legend("topright", inset =c(0.04,0.08), box.lty=0,xjust=0,yjust=0, legend=c("Couvain femelle", "Couvain m??le", "Taille population", "Energie miel"),col=c(" dark green", "dark blue","blue", "green"), lty=1,lwd=3, cex=0.7, title="Traits d??mographiques s??par??ment")
legend("topright", box.lty=0,inset =c(0.01,0.25),xjust=0,yjust=0,legend=c("Couvain*Energie miel", "Couvain femelle*Couvain m??le"), col=c("purple","black"), lty=5,lwd=3, cex=0.7, title="Traits d??mographiques en interaction")
legend("topright",box.lty=0,inset =c(0.07,0.35),xjust=0,yjust=0,legend=c("AFF", "LSP", "Activit?? de vols"), col=c("yellow","orange", "pink") ,lty=4,lwd=3, cex=0.7, title ="Traits de vie individuels")



