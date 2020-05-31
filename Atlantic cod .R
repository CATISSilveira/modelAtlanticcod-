# Alteracoes climaticas no bacalhau do atlantico (gadus morhua)
# Catarina Silveira a54337
# 31/05/2020
# Modelacao ecologica marinha e alteracoes climaticas
#https://github.com/CATISSilveira/modelAtlanticcod-

#setwd("C:/Users/Silveira/Desktop/courseMarineEcologicalModelling-master/Sessions [P]")

# Packages que e necessario usar para este script
#library(sp)
#library(rgdal)
#library(rgeos)
#library(raster)
#library(geosphere)
#library(ggplot2)
#library(gridExtra)
#library(rnaturalearth)
#library(rnaturalearthdata)
#library(leaflet)
#library(sf)
#library(sdmpredictors)
#source("sourceFunctions.R")


# Onde obti os dados 
world <- ne_countries(scale = 'medium')
recordsGBIF <- getOccurrencesGBIF("Gadus morhua")
recordsObis <- getOccurrencesObis("Gadus morhua")

records <- rbind(recordsGBIF,recordsObis)

# Remover as coordenadas que n?o tem l? nada
records <- removeNA(records,"Lon","Lat")

# Remover as coordenadas que est?o em duplicado 
records <- removeDuplicated(records,"Lon","Lat")

# Grafico depois de remover essas coordenadas
plot(world, col="Gray", border="Gray", axes=TRUE, main="Species distribution (gadus morhua)" , ylab="latitude", xlab="longitude")
points(records[,c("Lon","Lat")], pch=20, col="Black")
defineRegion(records,"Lon","Lat")
records <- selectRecords(records,"Lon","Lat")

# reconfirm that all records bellong to the know distribution of the species
plot(world, col="Gray", border="Gray", axes=TRUE, main="Species distribution (gadus morhua)" , ylab="latitude", xlab="longitude")
points(records[,c("Lon","Lat")], pch=20, col="Black")

# identify and remove records on land
records <- removeOverLand(records,"Lon","Lat")

# plot records V1
plot(world, col="Gray", border="Gray", axes=TRUE, main="Species distribution (gadus morhua)" , ylab="latitude", xlab="longitude")
points(records[,c("Lon","Lat")], pch=20, col="Black")

# plot records V2
ggplot() + 
  geom_polygon(data = world, aes(x = long, y = lat, group = group)) +
  geom_point(data = records, aes(x = Lon, y = Lat), color = "red") +
  scale_y_continuous(breaks = seq(-90,90, by=20)) +
  scale_x_continuous(breaks = seq(-180,180,by=20)) +
  coord_fixed()

# plot records V3
ggplot() + theme(
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank()) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group)) +
  geom_point(data = records, aes(x = Lon, y = Lat), color = "red") +
  scale_y_continuous(breaks = seq(-90,90, by=20)) +
  scale_x_continuous(breaks = seq(-180,180,by=20)) +
  coord_map("orthographic", orientation=c(20, 20, 0))

# save data frame to external file
Tabela1=write.table(records,file="",sep=";")

#..........................."........................................."............................"

# Segunda parte do trabalho predict para o presente 

layernames <- list_layers(datasets = "Bio-ORACLE")
View(layernames)

#Camadas que vou usar conforme a minha especie
environmentalConditions <- load_layers(c("BO2_tempmax_bdmean","BO2_tempmin_bdmean","BO2_dissoxmean_bdmean","BO2_ppmean_bdmean","BO2_salinitymean_bdmax"))
toKeep <- extract(environmentalConditions,records)
toKeep <- which(!is.na(toKeep[,1]))
records <- records[toKeep,]

# Coordenadas do atlantico uma vez que a minha especie pertence ao oceano atlantico 
europeanExtent <- extent(extent(-100,45,30.75,72.5))
environmentalConditions <- crop(environmentalConditions,europeanExtent)
background <- backgroundInformation(environmentalConditions,n=10000)
modelData <- prepareModelData(records,background,environmentalConditions)

folds <- get.block(records, background)
model <- train("Maxnet", modelData, folds = folds , fc="t" )

# given a set of possible hyperparameter values
h <- list(reg = seq(1,10,1) )
exp1 <- gridSearch(model, hypers = h, metric = "auc")
plot(exp1)
exp1@results

# fit a Maxent model to the dataset with the best hyperparameter values (o melhor valor 6)
model <- train("Maxnet", modelData, folds = folds , fc="t" , reg= 6)
getAUC(model, test = TRUE)

# determine relative variable contribution
viModel <- varImp(model, permut = 5)
viModel

# determine performance as AUC
getAUC(model, test = TRUE)

# reduce model complexity by dropping variabale at a time
reducedModel <- reduceVar(model, th = 5, metric = "auc", permut = 5)

# determine the final performance as AUC
getAUC(reducedModel, test = TRUE)

# determine the final relative variable contribution
viModel <- varImp(reducedModel, permut = 5)
viModel
plotVarImp(viModel)
# inspect response curves (Nao tenho o plot da producao primaria pois nao aparece no modelo reduzido)
plotResponse(reducedModel, var ="BO2_tempmax_bdmean", type = "logistic", only_presence = FALSE, marginal = FALSE, rug = FALSE, color="Black")
plotResponse(reducedModel, var ="BO2_dissoxmean_bdmean" , type = "logistic", only_presence = FALSE, marginal = FALSE, rug = FALSE, color="Black")
plotResponse(reducedModel, var ="BO2_salinitymean_bdmax" , type = "logistic", only_presence = FALSE, marginal = FALSE, rug = FALSE, color="Black")
plotResponse(reducedModel, var ="BO2_tempmin_bdmean" , type = "logistic", only_presence = FALSE, marginal = FALSE, rug = FALSE, color="Black")

# predict with Maxent to raster stack
mapPresent <- predict(reducedModel, environmentalConditions, type=c("logistic"))
plot(mapPresent)

# determine the threshold maximazing the sum of sensitivity and specificity
threshold <- thresholdMaxTSS(reducedModel)
threshold
getAccuracy(reducedModel,threshold = threshold)
getAUC(reducedModel, test = TRUE)

# generate a reclassification table
thresholdConditions <- data.frame(from = c(0,threshold) , to=c(threshold,1) , reclassValue=c(0,1))
thresholdConditions

# apply threshold to reclassify the predictive surface
mapPresentReclass <- reclassify(mapPresent, rcl = thresholdConditions)
plot(mapPresentReclass)

#.....................................".....................................".............................

#Prever para o futuro

# load layers of sea surface temperatures for the future
environmentalConditionsRCP26 <- load_layers(c("BO2_RCP26_2100_tempmax_bdmean","BO2_RCP26_2100_tempmin_bdmean","BO2_RCP26_2100_salinitymean_bdmax"))

#Esta layer nao se encontra nos packages 
rasterAdicional=raster("Data/rasterLayers/DissolvedMolecularOxygen Benthic Mean Pred Mean.tif")
environmentalConditionsRCP26 = stack(environmentalConditionsRCP26,rasterAdicional)

# change layer names for them to match
names(environmentalConditionsRCP26) <- c("BO2_tempmin_bdmean","BO2_tempmax_bdmean","BO2_dissoxmean_bdmean","BO2_salinitymean_bdmax")

# crop layers to the European extent
environmentalConditionsRCP26 <- crop(environmentalConditionsRCP26,europeanExtent)

# predict with Maxent to raster stack
mapRCP26 <- predict(reducedModel, environmentalConditionsRCP26, type=c("logistic"))
plot(mapRCP26)

