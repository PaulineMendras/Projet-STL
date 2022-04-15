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
plot.ts(ipi, xlab="Années", ylab="IPI", main="Production industrielle mensuelle de la fabrication
        de biscuits, biscottes et pâtisseries de conservation")
#La série ne parait pas stationnaire
#La série semble présenter une tendance descendante

#Autocorrélogramme total de la série
acf(ipi, lag.max=60)
#L'autocorrélation de premier ordre est très élevée et proche de 1 et l'autocorrélation totale décroit très lentement
#Cela confirme que la série n'est pas stationnaire

#Pour en être bien sûr, nous allons effectuer des tests de stationnarité
#Déterminons avant si la série présente une tendance déterministe, afin de choisir la spécification la mieux adaptée pour chaque test

trend <- 1:length(ipi)
lt <- lm(ipi ~ trend)
summary(lt)
#Le coefficient associé à la tendance est significatif au seuil de 1%
#La série comporte donc une tendance linéaire
#Le coefficient associé à la constante est également significatif au seuil de 1%

#Test de Dickey-Fuller augmenté, avec type constante + trend
adfTest(ipi, lag=0, type="ct")
#stationnaire...

#PP
summary(ur.pp(ipi, type="Z-tau", model="trend"))
#stationnaire...

#KPSS
summary(ur.kpss(ipi, type="tau", lags="long"))
#stationnaire...
