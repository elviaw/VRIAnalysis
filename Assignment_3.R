### GEOG 418/518 ###
### LAB 3 - SPATIAL AUTOCORRELATION ###

#If this is the first time you are running this code, you may have to install the following these packages that will be used this analysis. 
#To install them, please run the following script:  
install.packages("spdep")
install.packages("raster")
install.packages("rgdal")
install.packages("tmap")
install.packages("shinyjs")

#Load packages
library(rgdal)
library(tmap)
library(spdep)
library(raster)
library(shinyjs)

#Set your working folder directory, this folder will contain all data that R will read and output (graphs and tables, etc) will be written this is folder)
dir <- "/Users/epi/Documents/UVic/Geog518/Assignment 3/Data"
setwd(dir)

# Read data - Shapefile (do not include .shp)
VRI <- readOGR(dsn = ".", layer =  "WatershedVRI")
# Transform spatial reference of the original shapefile to NAD83 UTM Zone 10N
VRI <- spTransform(VRI, CRS("+init=epsg:26910")) #Find EPSG for NAD83 UTM Zone 10N

#view the first 6 rows of data
head(VRI@data)

#create a subset of the data that contain only the fields/columns that we are interested in 
vriCleanCols <- c("FID_VEG_CO", "POLYGON_ID", "PROJ_AGE_1",
                  "SITE_INDEX", "SPECIES__4", "SPECIES__5",
                  "PROJ_HEI_1", "SPECIES_PC", "SPECIES__6",
                  "VRI_LIVE_S", "BASAL_AREA", "WHOLE_STEM",
                  "CROWN_CL_1")
vriClean <- VRI[,vriCleanCols]

# Rename fields using  
newNames <- c("FID", "PolyID", "Stand_Age", "Site_Index",
              "CoDom_Sp", "Dom_Sp", "Stand_HT", "DomSP_Perc", 
              "CDomSP_Perc", "Stand_Dens", "Stand_BA", "Stand_StemBio", "Stand_CrownCl")
colnames(vriClean@data) <- newNames
#view first 6 rows of the data 
head(vriClean@data)

#Extract the attribute table of the shapefile as a data frame table 
vriCleanTable <- vriClean@data 

#Choose a variable of interest 
vriClean <- vriClean[!is.na(vriClean@data$Stand_HT), ]

#### We have cleaned the data and choose a variable to inspect (in this example, Height)

###Explore Data

# Explore map pallettes tool in R 
tmaptools::palette_explorer() #Note: Remember to close the palette_explorer window, R will not response to other command until this window is closed. 

# Visuallize the data using the pallettes tool 
map_StdHT <- tm_shape(vriClean) + 
  tm_polygons(col = "Stand_HT", 
              title = "Forest Height", 
              style = "jenks", 
              palette = "Greens", n = 7)

#view data in a map and in an interactive map
map_StdHT
tmap_mode("view")
tmap_mode("plot")

# Save the map in a .png format 
tmap_save(map_StdHT, "Forest Height.png", height=7)

# save the interactive map in html file
tmap_save(map_StdHT, "Forest Height.html")

### We have completed the section on visualizing/exploring the data 

### Assigning Neighborhood 

#Construct Queen's case neighborhood 
vri.nb <- poly2nb(vriClean)
#Vector files that will be used to assign weight
vri.net <- nb2lines(vri.nb, coords=coordinates(vriClean))
# Apply the same CRI to the data
crs(vri.net) <- crs(vriClean)

#View map 
Queen_Neighbour <- tm_shape(vriClean) + tm_borders(col='lightgrey') + 
  tm_shape(vri.net) + tm_lines(col='red')
# Save the map  
tmap_save(Queen_Neighbour, "Forest Height Neighborhood Queen.png", height=7)
tmap_save(Queen_Neighbour, "Forest Height Neighborhood Queen.html")

#Defining Neighbourhood using Rook's case
vri.nb2 <- poly2nb(vriClean, queen = FALSE)
vri.net2 <- nb2lines(vri.nb2, coords=coordinates(vriClean))
crs(vri.net2) <- crs(vriClean)

# View neighborhoods assigned by both Queen's case and Rook's case map in R (Plot and Viewer pane)
Neighbours <- tm_shape(vriClean) + tm_borders(col='lightgrey') + 
  tm_shape(vri.net) + tm_lines(col='blue', lwd = 2) +
  tm_shape(vri.net2) + tm_lines(col='yellow', lwd = 2) + 
  tm_add_legend('symbol', 
                col = c(adjustcolor( "blue", alpha.f = 1), adjustcolor( "yellow", alpha.f = 1)), shape = c(95, 95), 
                border.col = "grey40",
                labels = c('Queens case','Rooks case'),
                title="Neighbourhood Type")

# Save maps  
tmap_save(Neighbours, "Forest Height Neighbourhood Queen and Rook.png", height=7)
tmap_save(Neighbours, "Forest Height Neighbourhood Queen and Rook.html")

#Assign Spatial Weights (value) to the neighborhood, using style W, R will assigned value of 0 to 1 to the neighbours depending on the number of neighbours each point has, 
#Queen's Case
vri.lw <- nb2listw(vri.nb, zero.policy = TRUE, style = "W")
print.listw(vri.lw, zero.policy = TRUE)
#Rook's case
vri.lw2 <- nb2listw(vri.nb2, zero.policy = TRUE, style = "W")
print.listw(vri.lw2, zero.policy = TRUE)

###Completed the Neighborhood Matrix 

###Global Moran's I (Queen's case)
mi <- moran.test(vriClean$Stand_HT, vri.lw, zero.policy = TRUE)
mi 

#Calcualte the spatial range of the spatial autocorrelation 
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(vri.lw)

mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

# Calculate the z-value 
z <- (mI-eI)/sqrt(var) 
z

###Global Moran's I (Rook's case)
mi2 <- moran.test(vriClean$Stand_HT, vri.lw2, zero.policy = TRUE)
mi2 

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(vri.lw2)

#Assign values to calculate z value 
mI2 <- mi$estimate[[1]]
eI2 <- mi$estimate[[2]]
var2 <- mi$estimate[[3]]

z2 <- (1-eI)/sqroot(var) 
  
####Local Moran's I Analysis 
# Local Moran's I Analysis (Queen's Neighborhood), repeat for Rook's case
lisa.test <- localmoran(vriClean$Stand_HT, vri.lw, zero.policy = TRUE)
summary(lisa.test)

#Assign LISA result values to fields for mapping LIZA Z value 
vriClean$Ii <- lisa.test[,1]
vriClean$E.Ii<- lisa.test[,2]
vriClean$Var.Ii<- lisa.test[,3]
vriClean$Z.Ii<- lisa.test[,4]
vriClean$P<- lisa.test[,5]

#Mapping the LISA Z value results
map_LISAzvalue = subset(vriClean, P < 0.005) #Create new object with only P < 0.005 
map_LISAzvalue <- tm_shape(vriClean) +   tm_polygons(col = "lightgrey") +
  tm_shape(map_LISAzvalue) + 
  tm_polygons(col = "Z.Ii", 
              title = "Forest Height LISA Z Value (Queen's)", 
              style = "fisher", 
              palette = "RdBu", n = 9,) +
  tm_add_legend(type = "symbol", labels = c("Not significant"), col = c(adjustcolor("lightgrey", alpha.f = 1)), shape = c(22))
map_LISAzvalue
# Save the map  
tmap_save(map_LISAzvalue, "Forest Height LISA Z Value Queens.png", height=7)
tmap_save(map_LISAzvalue, "Forest Height LISA Z Value Queens.html")

#Plot the results of the LISA to a Scatter plot and save it as a png file 
png("LISA Moran Plot Queen's Case Height2.png")
mp <- moran.plot(vriClean$Stand_HT, vri.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Forest Height (location i) ", 
           ylab="Spatially Lagged Forest Height (neighbourhood of location i)", quiet=NULL)
if (require(ggplot2, quietly=TRUE)) {
  xname <- attr(mp, "xname")
  ggplot(mp, aes(x=x, y=wx)) + geom_point(shape=1) + 
    geom_smooth(formula=y ~ x, method="lm") + 
    geom_hline(yintercept=mean(mp$wx), lty=2) + 
    geom_vline(xintercept=mean(mp$x), lty=2) + theme_minimal() + 
    geom_point(data=mp[mp$is_inf,], aes(x=x, y=wx), shape=9) +
    geom_text(data=mp[mp$is_inf,], aes(x=x, y=wx, label=labels, vjust=1.5)) +
    xlab("Forest Height (location i)") + ylab("Spatially Lagged Forest Height (neighbourhood of location i)")
}
dev.off()

# Local Moran's I Analysis (Rooks's Neighborhood)
lisa.test2 <- localmoran(vriClean$Stand_HT, vri.lw2, zero.policy = TRUE)
summary(lisa.test2)

#Assign LISA result values to fields for mapping LIZA Z value 
vriClean$Ii2 <- lisa.test2[,1]
vriClean$E.Ii2<- lisa.test2[,2]
vriClean$Var.Ii2<- lisa.test2[,3]
vriClean$Z.Ii2<- lisa.test2[,4]
vriClean$P2<- lisa.test2[,5]

#Mapping the LISA Z value results
LISAzvalue2 = subset(vriClean, P2 < 0.005) #Create new object with only P < 0.005 
map_LISAzvalue2 <- tm_shape(vriClean) +   tm_polygons(col = "lightgrey") +
  tm_shape(LISAzvalue2) + 
  tm_polygons(col = "Z.Ii2", 
              title = "Forest Height LISA Z Value (Rook's)", 
              style = "fisher", 
              palette = "RdBu", n = 9,) +
  tm_add_legend(type = "symbol", labels = c("Not significant"), col = c(adjustcolor("lightgrey", alpha.f = 1)), shape = c(22))

map_LISAzvalue2
# Save the map  
tmap_save(map_LISAzvalue2, "Forest Height LISA Z Value Rooks.png", height=7)
tmap_save(map_LISAzvalue2, "Forest Height LISA Z Value Rooks.html")

#Plot the results of the LISA to a Scatter plot and save it as a png file 
png("LISA Moran Plot Rook's Case Height.png")
mp2 <- moran.plot(vriClean$Stand_HT, vri.lw2, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Forest Height (location i) ", 
                 ylab="Spatially Lagged Forest Height (neighbourhood of location i)", quiet=NULL)
if (require(ggplot2, quietly=TRUE)) {
  xname <- attr(mp2, "xname")
  ggplot(mp2, aes(x=x, y=wx)) + geom_point(shape=1) + 
    geom_smooth(formula=y ~ x, method="lm") + 
    geom_hline(yintercept=mean(mp2$wx), lty=2) + 
    geom_vline(xintercept=mean(mp2$x), lty=2) + theme_minimal() + 
    geom_point(data=mp2[mp2$is_inf,], aes(x=x, y=wx), shape=9) +
    geom_text(data=mp2[mp2$is_inf,], aes(x=x, y=wx, label=labels, vjust=1.5)) +
    xlab("Forest Height (location i)") + ylab("Spatially Lagged Forest Height (neighbourhood of location i)")
}
dev.off()
########################