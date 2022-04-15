########################################################################
################ PROJET DE SERIES TEMPORELLES LINEAIRES ################
################    Emeline LE HIR et Pauline MENDRAS   ################
########################################################################


rm(list=ls())

#Téléchargement des packages zoo, tseries et fUnitRoots
require(zoo)
require(tseries)
require(fUnitRoots)
require(urca)
require(vars)
require(stargazer)

################################
###  PARTIE 1 : LES DONNEES  ###
################################

##############
# Question 1 #
##############

#Importation des données
path <- "C:/Users/Pauline/Documents/GitHub/Projet-STL"
setwd(path)
datafile <- "valeurs_mensuelles.csv"
data <- read.csv(datafile,sep=";")

#Nettoyage des données
data <- data[-c(1:3),-3]
data<- data[seq(dim(data)[1],1),]
colnames(data)=c("Date","Indice")
data$Indice <- as.numeric(data$Indice)

##############
# Question 2 #
##############

#Représentation graphique de la série
ipi <- ts(data$Indice,start=c(1990,1), end=c(2022, 2), frequency = 12)
plot.ts(ipi, xlab="Années", ylab="IPI", main="Production industrielle mensuelle de la préparation
        de jus de fruits et de légumes")
#La série ne parait pas stationnaire
#La série semble présenter une tendance ascendante

#Autocorrélogramme total de la série
acf(ipi, lag.max=60)
#L'autocorrélation de premier ordre est très élevée et proche de 1 et l'autocorrélation totale décroit très lentement
#Cela confirme que la série n'est pas stationnaire

#Avant d'effectuer des tests de stationnarité pour confirmer la non-stationnarité de la série,
# déterminons si la série présente bien une tendance déterministe linéaire
trend <- 1:length(ipi)
lt <- lm(ipi ~ trend)
summary(lt)
#Le coefficient associé à la tendance est positif et significatif au seuil de 1%
#Le coefficient associé à la constante est significatif
#La série comporte donc une tendance linéaire

#Test de Dickey-Fuller augmenté
adf <- adfTest(ipi, lags=0, type="ct") #test ADF dans le cas avec constante et tendance

#Avant d'interpréter le test, vérifions que les résidus du modèle de régression sont bien non
#autocorrélés, sans quoi le modèle ne serait pas valide

#Test de validité des modèles: testent d'autocorrélation des résidus (comme dans le TD4)
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

#Comme la série est mensuelle, testons l'autocorrélation des résidus jusqu'à l'ordre 24 (2 ans)
Qtests(adf@test$lm$residuals,24,length(adf@test$lm$coefficients))
#L'absence d'autocorrélation des résidus est rejetée au moins 1 fois (Q(4) à Q(24)), le test ADF
#avec aucun retard n'est donc pas valide. Ajoutons des retards jusqu'à ce que les résidus
#ne soient plus autocorrélés

#Test validité des tests ADF comme en TD5
adfTest_valid <- function(series,kmax,type){ #tests ADF jusqu'à des résidus non autocorrélés
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")}
    else cat("nope \n")
    k <- k + 1
  }
  return(adf)
}

adf <- adfTest_valid(ipi,24,"ct")
#Il a fallu considérer 21 retards au test ADF pour supprimer l'autocorrélation des résidus

#Affichage des résultats du test valide maintenant
adf
#La racine unitaire n'est pas rejetée au seuil de 95% pour la série en niveau.
#La série est donc au moins I(1)

#1ère méthode pour rendre la série stationnaire: on garde la série corrigée de sa tendance
ipi_r <- lt$residuals 
plot.ts(r)
acf(r)
#La série ne semble pas stationnaire, nous ne sélectionnons pas cette méthode

#2ème méthode: différenciation première
diff_ipi <- ipi - lag(ipi,-1)
plot.ts(diff_ipi)
acf(diff_ipi, lag.max=60)
#La série semble sationnaire, vérifions cela avec des tests de stationnarité

#Avant cela, déterminons si la série présente une tendance déterministe pour choisir les spécifications adaptées aux tests
trend <- 1:length(diff_ipi)
lt <- lm(diff_ipi ~ trend)
summary(lt)
#Le coefficient associé à la tendance n'est pas significatif au seuil de 5%. La série ne présente donc pas de tendance linéaire.
#Le coefficient associé à la constante n'est pas sifnificatif au seuil de 5%

#Effectuons donc le test ADF dans le cas sans tendance ni constante, 
#en vérifiant l'absence d'autocorrélation des résidus, comme précédemment
adf_diff <- adfTest_valid(diff_ipi,24, type="nc")
#Il a fallu considérer 23 retards au test ADF pour supprimer l'autocorrélation des résidus
#On rejette l'hypothèse nulle de racine unitaire au seuil de 1%.
#Affichage des résultats du test valide maintenant
adf_diff
#Le test rejette la racine unitaire au seuil de 1%
#La série différenciée est donc stationnaire

#Test Phillips-Perron
summary(ur.pp(diff_ipi, type="Z-tau", model='constant'))
#On rejette l'hypothèse nulle de racine unitaire au seuil de 1%.
#La série différenciée est donc stationnaire

##############
# Question 3 #
##############

par(mfrow=c(2,1))
plot.ts(ipi, xlab="Années", ylab="IPI", main="Série avant transformation"); 
  plot.ts(diff_ipi, xlab="Années", ylab="diff_IPI", main="Série différenciée")

#################################
###  PARTIE 2 : MODELES ARMA  ###
#################################

##############
# Question 4 #
##############

#La série différenciée d'ordre 1 est stationnaire, donc d=1
#Pour choisir le modèle ARIMA(p,1,q) le plus approprié pour la série ipi (ARIMA(p,0,q) pour la série diff_ipi), on regarde l'autocorrélogramme total et l'autocorrélogramme partiel
diff_ipi <- diff_ipi - mean(diff_ipi) #on centre la série différenciée

par(mfrow=c(2,1))
acf(diff_ipi,60) #on considère qmax=1
pacf(diff_ipi,60) #on considère pmax=4
#On augmentera si les modèles retenus sont autocorrélés

#Les modèles possibles sont donc tous les modèles ARIMA(p,1,q) tels que p<=4 et q<=1
# Pour choisir le meilleur modèle, il faut vérifier la validité et l'ajustement
#des modèles possibles
#Pour tester la validité du modèle, nous utilisons la fonction Qtest définie précédemment
#Pour tester l'ajustement du modèle, nous regardons la significativité des coefficients
#Nous définissons pour cela la fonction signif comme dans le TD4
signif <- function(estim){ 
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

#Fonction d'affichage des tests pour la sélection du modèle ARIMA, comme dans le TD4
arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("tests de nullité des coefficients:\n")
  print(adjust)
  cat("\n tests d’absence d’autocorrélation des résidus:\n")
  print(pvals)
}

#Fonction pour estimer un ARIMA et en vérifier l'ajustement et la validité (TD4)
modelchoice <- function(p,q,data=diff_ipi, k=24) {
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,
                                          "masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,
           "ok"=ok))
}

#Fonction pour estimer et vérifier tous les ARIMA(p,q) avec p<=pmax et q<q=max (TD4)
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row){
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,")\n"))
    modelchoice(p,q)
  }))
}

pmax <- 4; qmax <- 1
armamodels <- armamodelchoice(pmax,qmax) #estime tous les ARIMA
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),]
#modèles bien ajustés et valides
selec
#Aucun modèle n'est ajusté et validé pour pmax=4 et qmax=1
#En reconsidérant l'ACF et la PACF, on prend le second ordre tel que p(h)<0 et r(h)<0, ie pmax=9 et qmax=8

pmax <- 13; qmax <- 13
armamodels <- armamodelchoice(pmax,qmax) #estime tous les ARIMA
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),]
#modèles bien ajustés et valides
selec
#2 modèles sont bien spécifiés: ARIMA(12,0,4) et ARIMA(4,0,12) pour la série diff_ipi

#Sélection des modèles via les  critères d'information
arima1204 <- arima(diff_ipi,c(12,0,4))
arima4012 <- arima(diff_ipi,c(4,0,12))
AIC(arima1204);BIC(arima1204)
AIC(arima4012);BIC(arima4012)
#Le modèle ARIMA(4,0,12) minimise les 2 critères d'infomation: nous retenons ce modèle

#Estimation du modèle ARIMA(4,0,12) pour la série diff_ipi
stargazer(arima4012, type = 'text')
summary(arima4012$residuals)
sigmaEpsilon <- sqrt(var(arima4012$residuals))
