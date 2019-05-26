#################################################################
# May 16 2019
# This script performs a basic batch "sanity" check analysis 
# ampinzonv@unal.edu.co
#
#################################################################

library(sybil)
library(glpkAPI)
library(sybilSBML)
library(minval)

##############################################################
# VARIABLES DE AMBIENTE (MODIFICAR DE ACUERDO A SU SISTEMA)
##############################################################
#Archivo de salida
f="_sanitycheck.csv"
outFile =paste(Sys.Date(),f,sep="")


#Ruta a carpeta con modelos
#modelsPath="/home/andres/allVMH"
modelsPath="/home/andres/calidadModelosVMH/models"

##############################################################
# RUTINAS DEL SCRIPT. EN ADELANTE NO MODIFICAR.
##############################################################

#
# The following function checks for blocked reactions in a model
# Originally scripted by Daniel Osorio.
#
blockedRxns <- function(model){
  bloqueadas <- NULL
  message (model@react_num, " reactions to be optimized. Can take a while!")
  for (reaccion in 1:model@react_num) {
    
    
    
    # Reinicia el coeficiente de objetivo a 0
    model@obj_coef <- rep(0, model@react_num)
    # Establece la reacci??n como funci??n objetivo
    model@obj_coef[reaccion] <- 1
    FBA <- optimizeProb(model)
    # Extrae los id's de las reacciones activas
    bloqueadas <- unique(c(bloqueadas, model@react_id[as.vector(FBA@fluxdist@fluxes!=0)]))
  }
  
  bloqueadas <- model@react_id[!model@react_id%in%bloqueadas]
  #just some feedback
  message (length(bloqueadas)," blocked reactions found!")
  
  return(length(bloqueadas))
  
}



####
#Inicializar el archivo de salida
####
header <- c("Name","Rxns","Mets","OptExitCode","SolStatus","objValue","OF_name","Genes","closedOptExitCode","closedSolStatus","closedOpjValue","blocked_rxns")
write(header,outFile, ncolumns=12, sep="\t",append = TRUE)

#Inicia iteración a tavés de todos los modelos en la carpeta definida
allModels <- list.files(path=modelsPath, pattern="*.xml", full.names=TRUE, recursive=FALSE)
for(i in 1:length(allModels)){
  
  model <- readSBMLmod(allModels[i])
  message ("=== Evaluating ",model@mod_name, " genome scale model ==")
  
  #Optimizar
  fbaSol <- optimizeProb(model, algorithm = "fba", retOptSol = TRUE)
  
 
  #Initialize the vector:populate it  with "ceros"
  results<-rep(0,12)
  #Guardar variables de optimización
  results[1]<-model@mod_name #Name
  results[2]<-model@react_num # Number of reactions in model
  results[3]<-model@met_num #Number of metabolites in model
  results[4]<-lp_ok(fbaSol) #exit code of optimization
  results[5]<-lp_stat(fbaSol) #Solution status of optimization
  results[6]<-lp_obj(fbaSol) #Optimal value of the objective function after optimization
  results[7]<-obj_func(fbaSol) #Objective function name
  results[8]<-length(model@genes) #Numero de genes en el modelo

 
 
 #Crear un modelo cerrado. Esto es que todas las reacciones de uptake son
 # puestas en 0. Este modelo se usará para ver si este es capaz de crecer
 # o generar productos aun sin fuentes de energia.
 exRxns <- findExchReact(model)
 upRxns <- uptReact(exRxns)
 closedModel <-changeBounds(model,upRxns,lb=0)
 
 #Optimizar nuevamente
 message("Closing the model and optimizing again.")
 closedFbaSol <- optimizeProb(closedModel, algorithm = "fba", retOptSol = TRUE)
  
 #Guardar resultados de la optimizacion del modelo cerrado
  results[9]<-lp_ok(closedFbaSol)
  results[10]<-lp_stat(closedFbaSol)
  results[11]<-lp_obj(closedFbaSol)
  results[12]<-blockedRxns(model)

##
## Save results. Write file.
##
write(results,outFile, ncolumns=12, sep="\t",append = TRUE)

#Unset variables, clean environment
#rm(list=ls())

} #ENDS ITERATION THROUGH ALL MODELS

