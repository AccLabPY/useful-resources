## Clase sobre aleatorización y poder estadístico
## Antonella Bandiera
## antonella.bandiera@itam.mx
##
## Realizado para la Red de Aprendizaje del AccLabPY
## 2021-12-16

rm(list=ls())

setwd("~/Dropbox/randomization_Paraguay_PNUD/clases")
library(tidyverse)
library(blockTools)
library(lmtest)
library(sandwich)


## Abrir los datos que vamos a usar para bloquear
covdata <- read.csv("final_covdata.csv")


dput(names(covdata))

## Elegir las variables que vamos a usar para bloquear

blockvars <- c("consejo_salud","sub_consejo_salud","rural","indigenas","poblacion","devf1", "devf2", 
               "devf3", "devf4", "healthf1", "healthf2","total_consultas", "not_same_block")
# armar los grupos: vamos a tener dos grupos estrictos (tesaireka) y dentro vamos a bloquear en 
# funciÃ³n de las otras

blockOut <- block(covdata,groups="tesaireka", #datos, grupos
                  n.tr = 2, # numero de tratamientos
                  id.vars="establecimientos", # el id 
                  block.vars = blockvars, # el vector con todas las variables que usamos para bloquear
                  algorithm = "optGreedy", # el algoritmo que arma los bloques
                  distance="mahalanobis") # la distancia a minimizar

# tenemos los dos grupos estrictos basados en el valor de la variable tesaireka
blockMat <- blockOut$blocks[[1]] 
blockMat2 <- blockOut$blocks[[2]]
# unimos a los dos grupos
blockMat <- rbind(blockMat, blockMat2)
covdata$blocknum <- NA
for(i in 1:nrow(blockMat)){
  blocknumUP <- as.numeric(rownames(blockMat[i,]))
  covdata$blocknum[as.character(covdata$establecimientos)%in%blockMat[i,c(1,2)]] <- blocknumUP
}

write.csv(covdata, 
          file="all_blocks2.csv",
          row.names=F)


# Vamos a generar la asignaciÃ³n, usando el enfoque re-randomization
# Idelamente queremos miles, para no demorar vamos a poner 100
#n.assign <- 2000
n.assign <- 100
assignments <- matrix(ncol=length(covdata$establecimientos)+2, nrow=n.assign)
colnames(assignments) <- c(as.character(covdata$establecimientos), "p","iter")
for(j in 1:n.assign){
  assignOut <- assignment(blockMat, seed=j)
  assignMat <- assignOut$assg[[1]]
  assignUp <- rep(0, nrow(covdata))
  for(i in 1:nrow(assignMat)){
    assignUp[as.character(covdata$establecimientos)%in%as.character(assignMat[i,1])] <- 1
  }
  fitUp <- lm(assignUp~as.matrix(covdata[,blockvars]))
  fitNull <- lm(assignUp~1)
  testpUp <- waldtest(fitUp, fitNull, vcov= vcovHC, test = c("F", "Chisq"))[[4]][2]
  assignments[j,1:length(covdata$establecimientos)] <- assignUp
  assignments[j,"iter"] <- j
  assignments[j,"p"] <- testpUp
  write.csv(assignments, file="assignments2.csv", row.names=F)
  cat(paste(j," ",sep="")); flush.console()
}

plot(ecdf(assignments[,"p"]))

#####
# Select from the restricted randomization 
######

assignmentIn <- read.csv("assignments2.csv")
plot(ecdf(assignmentIn$p))
abline(h=.5)
USFNames <- setdiff(names(assignmentIn), c("p","iter"))

#priorityUSF <- c("USF. UNION II - ESTANDAR", "USF.UNION I - ESTANDAR", "USF. CAPIIBARY II - ESTANDAR",
#                 "USF- CAPIIBARY I - ESTANDAR","USF- VILLA YGATYMI 1","USF- VILLA YGATIMI 2",
#                 "USF- YASY CANY I","USF- YASY CANY 2")
#assignmentIn$priorityProp <- apply(assignmentIn[,],
#                                 1,
#                                 function(x){mean(x[names(assignmentIn)%in%priorityUSF])}
#)
assignKeep <- subset(assignmentIn, p>.75#&priorityProp==.5
)
round(apply(assignKeep, 2, mean),2)
set.seed(121)
pickOne <- sample(1:nrow(assignKeep), 1)
pickOut <- t(as.matrix(assignKeep[pickOne,]))
write.csv(pickOut,
          file="pickOut2.csv",
          row.names=F)

## The one picked is n 43

assignments <- as.data.frame(assignments)
final_assignment <- assignments[which(assignments$iter==43), ]
final_assignment <- final_assignment[, -c(55:56)] # delete iter and p value

treatment_vector <- as.numeric(final_assignment[1,])
treatment_vector <- as.vector(treatment_vector)
covdata$treated <- NA
covdata$treated <- treatment_vector

write.csv(covdata, 
          file="covs_treatment_final.csv",
          row.names=F)



