## GEOMORPH EXAMPLE CODE ##

#LOAD LIBRARIES ----
#always do this first!!
library(geomorph) 
library(geiger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(gginnards)
library(ggphylomorpho)

#DATA IMPORT AND PREP ----

#Import data into R, make sure only first columns are not numbers
raw_data <- read_csv("Data/16minkeonly.csv") 
glimpse(raw_data)

#Average landmark takes
avg_takes <- raw_data %>% group_by(specimenID) %>% summarize(across(starts_with("Raw"), list(mean)))
glimpse(avg_takes)

#Create data frame with only numerical values and set row names as specimen ID
avg_takes <- remove_rownames(avg_takes)
has_rownames(avg_takes) #steps for tibbles
avg_takes <- data.frame(column_to_rownames(avg_takes, var = "specimenID")) #adding row names to data frame
glimpse(avg_takes)

#Save specimens names as object
specimens <- row.names(avg_takes) 

#Transform in 3D array, first number is number of landmarks, second is dimensions (2 or 3)
shape_array <- arrayspecs(avg_takes, 16, 3) 

##Extract classifier columns from raw data
#Make tibble for plots in ggplot
classifiers <- raw_data %>% 
  group_by(specimenID) %>% 
  summarize(age) %>%       #this creates a tibble with only the individual names and one text factor (e.g. age)
  distinct()               #this only keeps in the tibble unique combinations of the two values, so you get a tibble with the same number or rows as the data frame
glimpse(classifiers)
has_rownames(classifiers)    #always check

#Save classifiers as factors from tibble
factor_age <- as.factor(classifiers$age)

#remove objects that are not needed (e.g. "links")
remove(A) 


#GPA ALIGNMENT ----

#Procrustes alignment, should also show mean config coordinates
gpa<-gpagen(shape_array) 
plot(gpa) #see points in space

#Save Centroid size as object
Csize <- gpa$Csize 
#Log-transform Centroid size as object
logCsize <- log10(Csize) 

#Save mean shape to create links
mean_shape <- gpa$consensus 

#Coordinates of all specimens after GPA alignment
coords <- gpa$coords 

#Plot all specimens with mean to check that all landmarks are ok
plotAllSpecimens(coords, mean = TRUE, label = TRUE, plot.param = list(pt.cex = 0.5, mean.cex = 5, mean.bg = "black"))
#Save screenshot of 3D viewer
rgl.snapshot(filename = "Output/X.png") 

#Check for outliers, they would be displayed in red
plotOutliers(coords)

#Plot outliers shape against mean to check which landmarks ar the problem - use their number in the dataset OR save as object plotOutliers
plotRefToTarget(mean_shape,coords[,,4], method="TPS", mag = 1.5)
plotRefToTarget(mean_shape,coords[,,4], method="vector", mag = 1.5)


#PREPARE WARP MESH AND LINKS  ----

#Find specimen closer to mean, useful to create warp mesh
findMeanSpec(coords) #number below specimen name is the number of the specimen in the array

#Create object containing only that specimen coordinates
warp_specimen <- coords[,,2] #number displayed by findMeanSpec
warp_specimen 

#Import simplified mesh to create warp mesh on
mesh_3D <- read.ply("Data/simpleskull.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates to make sure it has same scale
range(mesh_3D$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(warp_specimen)

##Create warp mesh, to use as reference for visualization of analyses
ref_mesh <- warpRefMesh(mesh = mesh_3D, mesh.coord = warp_specimen, ref = mean_shape, color = NULL, centered = FALSE) 

##Define links to create skull shape - use 2D or easy configurations ONLY
links <- define.links(mean_shape, ptsize = 4) 
#check if table of links ok
View(links)
#to eliminate wrong links
links <- links[-24,] 
#see if links are ok in 3D space
plot(minke_gpa,links=links) 


#PCA COMPLETE DATASET ----

#Run PCA on complete dataset, color points by age/groups that are classifiers
PCA_all <- gm.prcomp(coords)

#List of PC components and proportion of variation
PCA_all 

##View plot
plot(PCA_all, main = "PCA all data",  pch = 21, #title and type of point to be used
     col = "deeppink",   #border of points
     bg = "deeppink",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_all$x[,1], y = PCA_all$x[,2], labels = rownames(PCA_all$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

#Save PC scores as object to use later
pcscores_all <- PCA_all$x 

#Save shapes of extremes for axes used in plot
PC1min_all <- PCA_all[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_all <- PCA_all[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_all <- PCA_all[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_all <- PCA_all[["shapes"]][["shapes.comp2"]][["max"]] 

#Show deformation grids on axis from mean shape, do this for all 4 extremes - "TPS" method
plotRefToTarget(mean_shape, PC1min_all, method = "TPS", mag = 1, label = FALSE)     #save image

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - "points" method
plotRefToTarget(mean_shape, PC1min_all, method = "points", mag = 1, label = FALSE)    #save as HTML

#Show 3D deformation from mean with vectors, do this for all 4 extremes - "vector" method
plotRefToTarget(mean_shape, PC1min_all, method = "vector", mag = 1, label = FALSE)    #save as screenshot

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(mean_shape, PC1min_all, mesh = ref_mesh, method = "surface", mag = 1, label = FALSE)   #save as HTML
 
##3D windows save
#Save screenshot of 3D window, useful for lateral and dorsal views - use screen snip if it fails
rgl.snapshot(filename = "Output/X.png") 
#Save 3D window as html file - 3D widget
scene <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_all <- as_tibble(pcscores_all)
glimpse(pcscores_all)
#Add labels and other attributes to tibble as columns
pcscores_all <- pcscores_all %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(pcscores_all)

##Order tibble by variable e.g. age
#Make factor for variable
pcscores_all$age <- factor(pcscores_all$age, 
                    levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
pcscores_all <- pcscores_all_size[order(pcscores_all$age),]
glimpse(pcscores_all)

#Nice plot
ggplot(pcscores_all, aes(x = Comp1, y = Comp2, label = individuals, colour = age))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+            #legend and color adjustments
  theme_bw()+
  xlab("PC 1 (48.89%)")+ #copy this from standard PCA plot
  ylab("PC 2 (22.17%)")+
  ggtitle("PCA all data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position

##Regression PC1 and PC2 vs logCS
#Create data frame with data
pcscores_all_size <- data.frame(PCA_all$x , size = logCsize)

#Calculate regression for each component
reg_PC1all_size <- lm(Comp1 ~ size, data = pcscores_all_size)
reg_PC2all_size <- lm(Comp2 ~ size, data = pcscores_all_size)

#View results and p-value
summary(reg_PC1all_size)
summary(reg_PC2all_size)

#Diagnostic plots for both regressions
init <- par(no.readonly=TRUE) #store initial plot parameters to restore later
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
plot(reg_PC1all_size, cex = 1.2, font.main = 2)
plot(reg_PC2all_size, cex = 1.2, font.main = 2)                             
par(init)                     #restore initial plot parameters (1 plot showing at a time)

#Plot regression with ggplot
#Convert data frame to tibble
pcscores_all_size <- as_tibble(pcscores_all_size)
glimpse(pcscores_all_size)
#Add labels and other attributes to tibble as columns
pcscores_all_size <- pcscores_all_size %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(pcscores_all_size)

#Order tibble by variable e.g. age
#Make factor for variable
pcscores_all_size$age <- factor(pcscores_all_size$age, 
                          levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
pcscores_all_size <- pcscores_all_size[order(pcscores_all_size$age),]
glimpse(pcscores_all_size)

#Create data frame with line parameters from regression
#Allows to show a line on PC plot with specimens colored IF col and other graphics OUTSIDE of aes()!!!
regline_PCAall_size <- data.frame(int1 = reg_PC1all_size[["coefficients"]][["(Intercept)"]], 
                                  slope1 = reg_PC1all_size[["coefficients"]][["size"]], #PC1 values
                                  int2 = reg_PC2all_size[["coefficients"]][["(Intercept)"]], 
                                  slope2 = reg_PC2all_size[["coefficients"]][["size"]]) #PC2 values
regline_PCAall_size 

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(pcscores_all_size, aes(x = size, y = Comp1, label = individuals, colour = age))+
  #line on plot
  geom_abline(data = regline_PCAall_size, aes(intercept = int1, slope = slope1), colour = "darkblue", size = 0.8, linetype = "dashed", show.legend = FALSE)+
  geom_point(size = 3)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+           
  theme_classic(base_size = 12)+
  xlab("logCS")+
  ylab("PC 1 (48.89%)")+
  ggtitle("PC1 vs logCS - p-value <0.05**")+  #copy from summary linear model
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))+
  geom_text_repel(colour = "black", size = 3.5)

#Repeat for other component
#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(pcscores_all_size, aes(x = size, y = Comp2, label = individuals, colour = age))+
  #line on plot
  geom_abline(data = regline_PCAall_size, aes(intercept = int2, slope = slope2), colour = "darkblue", size = 0.8, linetype = "dashed", show.legend = FALSE)+
  geom_point(size = 3)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+           
  theme_classic(base_size = 12)+
  xlab("logCS")+
  ylab("PC 2 (22.17%)")+
  ggtitle("PC2 vs logCS - p-value 0.2468")+  #copy from summary linear model
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))+
  geom_text_repel(colour = "black", size = 3.5)


#ALLOMETRY CORRECTION ----
##Evaluate allometry and get the allometry-free shapes using LogCS, use this for analyses

#Regression shape on logCS size
allometry <- procD.lm(coords ~ logCsize, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry) 

#Create residuals array to then save as coordinates for analyses
allometry_array <- arrayspecs(allometry$residuals,p = dim(coords)[1], k = dim(coords)[2]) 

#New shapes adjusted for allometry with CS to use in analyses
allometry_residuals <- allometry_array + array(mean_shape, dim(allometry_array)) 

#Save mean shape of allometry-adjusted shapes to sue later
mean_shape_residuals <- mshape(allometry_residuals)

##Plot shape vs logCS to visualize allometry
#Diagnostic plots to check if model is appropriate - similar to ANOVA tables
init <- par(no.readonly=TRUE) #store initial plot parameters to restore later
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
plot(allometry,type = "diagnostics", cex = 1.2, font.main = 2)
par(init)                     #restore initial plot parameters (1 plot showing at a time)

#Regression score of shape vs logCS - regression method with "RegScore" plotting
allometry_plot_regscore <- plot(allometry, type = "regression",predictor = logCsize, reg.type = "RegScore",
                          main = "Shape vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics
text(x = logCsize, y = allometry_plot_regscore$RegScore, labels = specimens,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

##Add regression line with confidence intervals to plot
#Create object to use for linear model
allometry_regscores <- allometry_plot_regscore[["RegScore"]] 

#Linear model for line
allometry_regline <- (lm(allometry_regscores ~ logCsize))

#Draw line on plot
abline(allometry_regline, col = "darkseagreen3", 
       lty = 2, lwd = 2)     #line type (e.g. dashed) and width

#Add confidence intervals
#Create data for confidence intervals
allometry_xvals <- seq(min(logCsize), max(logCsize), length = 12)   #use min and max of x values (logCS) as limits and use number of specimens as length of sequence
allometry_newX <- expand.grid(logCsize = allometry_xvals)                     #warp x_vals on values of x axis (logCS)
allometry_newY <- predict(allometry_regline, newdata = data.frame(x = allometry_newX), interval="confidence",
        level = 0.95)                                      #predict the y values based on the x sequence

#Draw confidence intervals lines on plot
matlines(allometry_newX, allometry_newY[,2:3],                                 #first column of newY not useful, it is the fit, 2 and 3 are the min and max values
         col = "darkseagreen3", lty=1)                     #line graphics

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_plot_tibble <- data.frame(logCS = allometry_plot_regscore[["plot.args"]][["x"]], RegScores = allometry_plot_regscore[["plot.args"]][["y"]])
allometry_plot_tibble

#Convert data frame to tibble
allometry_plot_tibble <- as_tibble(allometry_plot_tibble)
glimpse(allometry_plot_tibble)
#Add labels and other attributes to tibble as columns
allometry_plot_tibble <- allometry_plot_tibble %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(allometry_plot_tibble)

#Order tibble by variable e.g. age
#Make factor for variable
allometry_plot_tibble$age <- factor(allometry_plot_tibble$age, 
                                       levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
allometry_plot_tibble <- allometry_plot_tibble[order(allometry_plot_tibble$age),]
glimpse(allometry_plot_tibble)

##Add regression line with confidence intervals
#Make data frame of data for confidence intervals
allometry_conf_intervals <- data.frame(allometry_newX, allometry_newY)
allometry_conf_intervals 
#Rename columns to match main plot tibble variables for x and y
allometry_conf_intervals <- rename(allometry_conf_intervals, logCS = logCsize, RegScores = fit)
allometry_conf_intervals 

#Convert data frame to tibble
allometry_conf_intervals <- as_tibble(allometry_conf_intervals)
glimpse(allometry_conf_intervals)
#Add labels and other attributes to tibble as columns to match main plot tibble
allometry_conf_intervals <- allometry_conf_intervals %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(allometry_conf_intervals)

#Order tibble by variable e.g. age to match main plot tibble
#Make factor for variable
allometry_conf_intervals$age <- factor(allometry_conf_intervals$age, 
                                levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
allometry_conf_intervals <- allometry_conf_intervals[order(allometry_conf_intervals$age),]
glimpse(allometry_conf_intervals)

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(allometry_plot_tibble, aes(x = logCS, y = RegScores, label = individuals, colour = age))+
  geom_smooth(data = allometry_conf_intervals, aes(ymin = lwr, ymax = upr), stat = 'identity',          #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+           
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Shape vs logCS - p-value = 0.004**")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1)     #position of tables relative to point (proximity and distance) 


#PC1 values vs logCS - regression method with "PredLine" plotting
allometry_plot_predline <- plot(allometry,type = "regression",predictor = logCsize, reg.type = "PredLine",
                             main = "PC1 vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)
text(x = logCsize, y = allometry_plot_predline$PredLine, labels = specimens,
     pos = 2, offset = 0.5, cex = 0.75)  

#PCA plot of fitted values - likely not very useful if size is a big component of variation, only PC1 will have weight
allometry_plot_PC <- plot(allometry,type = "PC",
                         main = "PCA fitted values", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)
text(x = allometry_plot_PC$plot.args$x, y = allometry_plot_PC$plot.args$y, labels = specimens,
     pos = 3, offset = 0.5, cex = 0.75) 

#ALLOMETRY CAC PLOT ----
##Specific plotAllometry functions - "CAC" and "size.shape" PCA (unclear what is the purpose of size-shape PCA, will not use)

#CAC plot - for simple allometry is the same as RegScore plot and two-block PLS plot of shape and size
allometry_plot_CAC <- plotAllometry(allometry, size = logCsize, logsz = FALSE, method = "CAC", 
                      main = "CAC plot", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics

##Make better CAC allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_plot_CAC_tibble  <- data.frame(logCS = allometry_plot_CAC[["size.var"]], CAC = allometry_plot_CAC[["CAC"]], #common allometric component
                                      RSC1 = allometry_plot_CAC[["all.plot.args"]][["RSC"]][["x"]])               #main residual shape component
allometry_plot_CAC_tibble

#Convert data frame to tibble
allometry_plot_CAC_tibble <- as_tibble(allometry_plot_CAC_tibble)
glimpse(allometry_plot_CAC_tibble)
#Add labels and other attributes to tibble as columns
allometry_plot_CAC_tibble <- allometry_plot_CAC_tibble %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(allometry_plot_CAC_tibble)

#Order tibble by variable e.g. age
#Make factor for variable
allometry_plot_CAC_tibble$age <- factor(allometry_plot_CAC_tibble$age, 
                                       levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
allometry_plot_CAC_tibble <- allometry_plot_CAC_tibble[order(allometry_plot_CAC_tibble$age),]
glimpse(allometry_plot_CAC_tibble)

##Two plots: CAC vs logCsize (A) and RSC1 vs CAC (B)
#Plot A
##Add regression line with confidence intervals
#Create object to use for linear model
allometry_CAC_regscores <- allometry_plot_CAC[["CAC"]]

#Linear model for line
allometry_CAC_regline <- (lm(allometry_CAC_regscores ~ logCsize))

#Add confidence intervals
#Create data for confidence intervals
allometry_CAC_xvals <- seq(min(logCsize), max(logCsize), length = 12)   #use min and max of x values (logCS) as limits and use number of specimens as length of sequence
allometry_CAC_newX <- expand.grid(logCsize = allometry_CAC_xvals)                     #warp x_vals on values of x axis (logCS)
allometry_CAC_newY <- predict(allometry_CAC_regline, newdata = data.frame(x = allometry_CAC_newX), interval="confidence",
                level = 0.95)                                      #predict the y values based on the x sequence

#Make data frame of data for confidence intervals
allometry_CAC_conf_intervals <- data.frame(allometry_CAC_newX, allometry_CAC_newY)
#Rename columns to match main plot tibble varibales for x and y
allometry_CAC_conf_intervals <- rename(allometry_CAC_conf_intervals, logCS = logCsize, CAC = fit)
allometry_CAC_conf_intervals

#Convert data frame to tibble
allometry_CAC_conf_intervals <- as_tibble(allometry_CAC_conf_intervals)
glimpse(allometry_CAC_conf_intervals)
#Add labels and other attributes to tibble as columns to match main plot tibble
allometry_CAC_conf_intervals <- allometry_CAC_conf_intervals %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(allometry_CAC_conf_intervals)

#Order tibble by variable e.g. age
#Make factor for variable
allometry_CAC_conf_intervals$age <- factor(allometry_CAC_conf_intervals$age, 
                                levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
allometry_CAC_conf_intervals <- allometry_CAC_conf_intervals[order(allometry_CAC_conf_intervals$age),]
glimpse(allometry_CAC_conf_intervals)

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(allometry_plot_CAC_tibble, aes(x = logCS, y = CAC, label = individuals, colour = age))+
  geom_smooth(data = allometry_CAC_conf_intervals, aes(ymin = lwr, ymax = upr), stat = 'identity',          #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+           
  theme_classic(base_size = 12)+
  ggtitle ("CAC vs logCS - p-value = 0.004**")+ #copy from allometry model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
               force_pull = 3, point.padding = 1) 

#Plot B

#Nice plot with specimens colored by age - no labels
ggplot(allometry_plot_CAC_tibble, aes(x = CAC, y = RSC1, colour = age))+
 #confidence intervals and reg line using standard function, difficult to do with external model for residuals (too small)
 #see function in CAC plot A for code for external model
  geom_smooth(method = 'lm',        #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.5)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+           
  theme_classic(base_size = 12)+
  ggtitle ("Residual shape component (RSC) 1 vs CAC")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


#TWO-BLOCK PLS ----

#Two-block PLS between logCS and shape, another way to visualize allometry
shape_logCS_pls <- two.b.pls(logCsize, coords,iter = 999) 

#Get P-value of regression
shape_logCS_pls

#Plot two-block PLS with regression line
shape_logCS_pls_plot <- plot(shape_logCS_pls, 
      pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2)   #improve appearance of points
#Save plot arguments as objects to use in plots
block1_pls_logCS <- shape_logCS_pls_plot$plot.args$x
block2_pls_shape <- shape_logCS_pls_plot$plot.args$y
#Add labels
text(x = block1_pls_logCS, y = block2_pls_shape, labels = specimens,
     pos = 3, offset = 0.5, cex = 0.75)   #improve appearance of labels

##Make better PLS plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
shape_logCS_pls_plot_tibble <- data.frame(block1 = block1_pls_logCS, block2 = block2_pls_shape)
shape_logCS_pls_plot_tibble

#Convert data frame to tibble
shape_logCS_pls_plot_tibble <- as_tibble(shape_logCS_pls_plot_tibble)
glimpse(shape_logCS_pls_plot_tibble)
#Add labels and other attributes to tibble as columns
shape_logCS_pls_plot_tibble <- shape_logCS_pls_plot_tibble %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(shape_logCS_pls_plot_tibble)

#Order tibble by variable e.g. age for plot legend
#Make factor for variable
shape_logCS_pls_plot_tibble$age <- factor(shape_logCS_pls_plot_tibble$age, 
                                levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
shape_logCS_pls_plot_tibble <- shape_logCS_pls_plot_tibble[order(shape_logCS_pls_plot_tibble$age),]
glimpse(shape_logCS_pls_plot_tibble)

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(shape_logCS_pls_plot_tibble, aes(x = block1, y = block2, label = individuals, colour = age))+
#confidence intervals and reg line, before points  
  geom_smooth(method='lm',     #use standard function, values too small    
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+  
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), 
                        values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+          
  theme_classic(base_size = 12)+
  xlab("PLS1 Block 1: logCS")+
  ylab("PLS1 Block 2: Shape")+
  ggtitle ("PLS logCS vs Shape - p-avlue = 0.008**")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5), axis.title = element_text(size = 11))


#PCA ALLOMETRY RESIDUALS ----

#New PCA plot with data corrected for allometry
PCA_residuals <- gm.prcomp(allometry_residuals) 

#List of PC components and proportion of variations
PCA_residuals

##View plot
plot(PCA_residuals, main = "PCA residuals",  pch = 21, #title and type of point to be used
     col = "deeppink", bg = "deeppink", cex = 1, font.main = 2)  #improve graphics
#Add quick labels to plot
text(x = PCA_residuals$x[,1], y = PCA_residuals$x[,2], labels = rownames(PCA_residuals$x), 
     pos = 1, offset = 0.5, cex = 0.75)    #improve graphics

#Save PC scores as object to use later
pcscores_res <- PCA_residuals$x

#Save shapes of extremes for axes used in plot
PC1min_res <- PCA_residuals[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_res <- PCA_residuals[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_res <- PCA_residuals[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_res <- PCA_residuals[["shapes"]][["shapes.comp2"]][["max"]] 

#Show deformation grids on axis from mean shape, do this for all 4 extremes - "TPS" method
plotRefToTarget(mean_shape_residuals, PC1min_res, method = "TPS", mag = 1, label = FALSE)  #save image

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - "points" method
plotRefToTarget(mean_shape_residuals, PC1min_res, method = "points", mag = 1, label = FALSE)   #save as HTML

#Show 3D deformation from mean with vectors, do this for all 4 extremes - "vector" method
plotRefToTarget(mean_shape_residuals, PC1min_res, method = "vector", mag = 1, label = FALSE)   #save as screenshot

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(mean_shape_residuals, PC1min_res, mesh = ref_mesh, method = "surface", mag = 1, label = FALSE)   #save as HTML

##3D windows save
#Save screenshot of 3D window, useful for lateral and dorsal views - use screen snip if it fails
rgl.snapshot(filename = "Output/X.png") 
#Save 3D window as html file - 3D widget
scene <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_res <- as_tibble(pcscores_res)
glimpse(pcscores_res)
#Add labels and other attributes to tibble as columns
pcscores_res <- pcscores_res %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(pcscores_res)

#Nice plot
ggplot(pcscores_res, aes(x = Comp1, y = Comp2, label = individuals, colour = age))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), 
                        values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+           #legend and color adjustments
  theme_bw()+
  xlab("PC 1 (44.15%)")+ #copy this from standard PCA plot
  ylab("PC 2 (22.79%)")+
  ggtitle("PCA residuals")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position


#ANOVA and MORPHOLOGICAL DISPARITY AMONG GROUPS  ----

#Create dataframe to operate more easily - use allometry residuals
allometry_df <- geomorph.data.frame(allometry_residuals, gp = factor_age) 

#Conduct ANOVA to test if there is significant shape variatiojn among groups - are shapes in each group different from the otehr groups?
group_anova <- procD.lm(allometry_residuals ~ gp, iter=999, data = allometry_df) 

#Results and significance of ANOVA
summary(group_anova) 

#Diagnostic plots to check if model is appropriate - similar to ANOVA tables - DO also if not significant, don't do other plots
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
plot(group_anova, type = "diagnostics", cex = 1.2, font.main = 2)
par(init)                     #restore initial plot parameters (1 plot showing at a time)


##Morphological disparity 
#Calculate Procrustes variances and distances between groups after allometric correction - how different are each group shapes compared to other group shapes?
group_disparity <- morphol.disparity(allometry_residuals ~ 1, groups = ~ factor_age, iter = 999) #do NOT use dataframe, issues with groups

#Results and significance
summary(group_disparity) 

#TRAJECTORY ANALYSIS ----
#Shows trajectories of variation using groups, use obj from procD.lm

#Trajectory points must be defined and they can either be a factor to assess different trajectories within groups or they should be = to the nummber of groups
group_trajectory <- trajectory.analysis(allometry, groups = classifiers$age, traj.pts = 4, pca = TRUE, print.progress = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different?
summary(group_trajectory, show.trajectories = TRUE, attribute = "MD") 
#Trajectory correlations - only useful if there are significant differences 
summary(group_trajectory, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
#Trajectory shape differences - can only be used if each trajectory/group has 3 or more points/specimens - what is the distance between trajectories?
summary(group_trajectory, show.trajectories = TRUE, attribute = "SD") 

#Plot results - PCA of fitted values, similar to CVA
group_trajectory_plot <- plot(group_trajectory, main = "Trajectories by age",  pch = 21, #title and type of point to be used
                        col = "gray", bg = "gray", cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(group_trajectory_plot, 
                 traj.pch = 21, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = 3, end.bg = 2) #trajectory line graphics

#doesn't work -> legend("right", legend = factor_age)

##Make better PCA plot using ggplot
#Read PC scores as tibble
group_trajectory_pcscores <- as_tibble(group_trajectory_plot[["pc.points"]])
glimpse(group_trajectory_pcscores)

#Find out order of group variables - should be ordered by factor as done when creating classifiers tibble at the start
glimpse(group_trajectory_plot[["trajectories"]])

#Create data frame that contains the group variables in order and that has = number of rows to pc score points
#Needed to find end points of trajectory lines for plot
rows_to_repeat <- bind_rows(data.frame(group = "earlyFetus", counter = 1:12), data.frame(group = "lateFetus", counter = 13:24), 
                            data.frame(group = "neonate", counter = 25:36), data.frame(group = "adult", counter = 37:48))
glimpse(rows_to_repeat)
#Delete extra counter column
rows_to_repeat <- rows_to_repeat [,-2]

#Add group names and other attributes to tibble as columns
group_trajectory_pcscores <- group_trajectory_pcscores %>% mutate(age = rows_to_repeat)
glimpse(group_trajectory_pcscores)

#Order tibble by variable e.g. age
#Make factor for variable
group_trajectory_pcscores$age <- factor(group_trajectory_pcscores$age, 
                                         levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
group_trajectory_pcscores <- group_trajectory_pcscores[order(group_trajectory_pcscores$age),]
glimpse(group_trajectory_pcscores)

#Calculate means of PC1 and PC2 per group to draw trajectories
group_trajectory_pcscores_means <- group_trajectory_pcscores %>% group_by(age)  %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(group_trajectory_pcscores_means)

#Order tibble by variable e.g. age
#Make factor for variable
group_trajectory_pcscores_means$age <- factor(group_trajectory_pcscores_means$age, 
                                       levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
group_trajectory_pcscores_means <- group_trajectory_pcscores_means[order(group_trajectory_pcscores_means$age),]
glimpse(group_trajectory_pcscores_means)


#Rename columns so they are easier to use for plot
group_trajectory_pcscores_means <- group_trajectory_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
group_trajectory_pcscores_means

#Nice plot
ggplot(group_trajectory_pcscores, aes(x = PC1, y = PC2, colour = age))+
  geom_point(size = 3)+
#add trajectory lines, one line for each, write row number from tibble, should be in order as legend of plot
     geom_segment(data = group_trajectory_pcscores_means, aes(x = x[1], y = y[1],  #earlyF
                                                       xend =  x[2], yend = y[2], colour = age), #lateF
                 colour = "snow4", size = 0.8, linejoin = 'mitre', arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))+ #add arrow at end  
     geom_segment(data = group_trajectory_pcscores_means, aes(x = x[2], y = y[2], #lateF
                                                       xend = x[3], yend = y[3], colour = age), #neonate
                 colour = "snow4", size = 0.8, linejoin = 'mitre', arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))+
     geom_segment(data = group_trajectory_pcscores_means, aes(x = x[3], y = y[3], #neonate
                                                       xend = x[4], yend = y[4]),  #adult 
                 colour = "snow4", size = 0.8, linejoin = 'mitre', arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+            #legend and color adjustments
  theme_bw()+
  xlab("PC 1 (60.17%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (30.92%)")+
  ggtitle("Trajectories by age")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position


#SYMMETRY ANALYSIS ----

#Create links between pairs of symmetric landmarks across line of symmetry to then run symmetry analysis
landpairs <- define.links(mean_shape, ptsize = 4) 
landpairs    #could be imported and then transformed in a 2D array

#Analysis of bilateral symmetry to then work on symmetric component alone
skull_symmetry <- bilat.symmetry(coords, ind = specimens, object.sym = TRUE, land.pairs = landpairs) 

#Check ANOVA results - significance value of analysis - if no significant value for "side" means objects are symmetrical
summary(skull_symmetry)

#Create object to work on symmetric component alone
symmetry_residuals <- skull_symmetry$symm.shape 

#Mean shape of symmetric component, to use for analysis
mean_shape_symmetry <- mshape(symmetry_residuals) 

#Plot on mesh to visualize asymmetry
plot(skull_symmetry, warpgrids = TRUE, mesh = ref_mesh) 
#Save 3D window as html file - 3D widget
scene <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

##Perform allometric correction and PCA on symmetry-only shapes
##Allometry regression
allometry_sym <- procD.lm(symmetry_residuals ~ logCsize, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry_sym) 

#Create residuals array to then save as coordinates for analyses
allometry_residuals_sym <- arrayspecs(allometry_sym$residuals, p = dim(symmetry_residuals)[1], k = dim(symmetry_residuals)[2]) 

#New shapes adjusted for allometry with CS to use in analyses
allometry_residuals_sym <- allometry_residuals_sym + array(mean_shape_symmetry, dim(allometry_residuals_sym)) 

#Save mean shape of allometry-adjusted shapes to sue later
mean_shape_residuals_sym <- mshape(allometry_residuals_sym)

##Plot shape vs logCS to visualize allometry
#Diagnostic plots to check if model is appropriate - similar to ANOVA tables
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
plot(allometry_sym, type = "diagnostics", cex = 1.2, font.main = 2)
par(init)                     #restore initial plot parameters (1 plot showing at a time)


#Regression score of shape vs logCS - regression method with regression score plotting
allometry_sym_plot_regscore <- plot(allometry_sym,type = "regression",predictor = logCsize, reg.type = "RegScore",
                               main = "Symmetric shape vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_sym_plot_tibble <- data.frame(logCS = allometry_sym_plot_regscore[["plot.args"]][["x"]], 
                                      RegScores = allometry_sym_plot_regscore[["plot.args"]][["y"]])

#Convert data frame to tibble
allometry_sym_plot_tibble <- as_tibble(allometry_sym_plot_tibble)
glimpse(allometry_sym_plot_tibble)
#Add labels and other attributes to tibble as columns
allometry_sym_plot_tibble <- allometry_sym_plot_tibble %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(allometry_sym_plot_tibble)

#Order tibble by variable e.g. age
#Make factor for variable
allometry_sym_plot_tibble$age <- factor(allometry_sym_plot_tibble$age, 
                                       levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
allometry_sym_plot_tibble <- allometry_sym_plot_tibble[order(allometry_sym_plot_tibble$age),]
glimpse(allometry_sym_plot_tibble)

##Add regression line with confidence intervals
#Create object to use for linear model
allometry_sym_regscores <- allometry_sym_plot_regscore[["RegScore"]] 

#Linear model for line
allometry_sym_regline <- (lm(allometry_sym_regscores ~ logCsize))

#Add confidence intervals
#Create data for confidence intervals
allometry_sym_xvals <- seq(min(logCsize), max(logCsize), length = 12)   #use min and max of x values (logCS) as limits and use number of specimens as length of sequence
allometry_sym_newX <- expand.grid(logCsize = allometry_sym_xvals)                     #warp x_vals on values of x axis (logCS)
allometry_sym_newY <- predict(allometry_sym_regline, newdata = data.frame(x = allometry_sym_newX), interval="confidence",
                level = 0.95)                                      #predict the y values based on the x sequence

#Make data frame of data for confidence intervals
allometry_sym_conf_intervals <- data.frame(allometry_sym_newX, allometry_sym_newY)
#Rename columns to match main plot tibble varibales for x and y
allometry_sym_conf_intervals <- rename(allometry_sym_conf_intervals, logCS = logCsize, RegScores = fit)
allometry_sym_conf_intervals

#Convert data frame to tibble
allometry_sym_conf_intervals <- as_tibble(allometry_sym_conf_intervals)
glimpse(allometry_sym_conf_intervals)
#Add labels and other attributes to tibble as columns to match main plot tibble
allometry_sym_conf_intervals <- allometry_sym_conf_intervals %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(allometry_sym_conf_intervals)

#Order tibble by variable e.g. age
#Make factor for variable
allometry_sym_conf_intervals$age <- factor(allometry_sym_conf_intervals$age, 
                                       levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
allometry_sym_conf_intervals <- allometry_sym_conf_intervals[order(allometry_sym_conf_intervals$age),]
glimpse(allometry_sym_conf_intervals)


#Nice plot with regression line and confidence intervals
ggplot(allometry_sym_plot_tibble, aes(x = logCS, y = RegScores, label = individuals, colour = age))+
  geom_smooth(data = allometry_sym_conf_intervals, aes(ymin = lwr, ymax = upr), stat = 'identity',          #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), 
                      values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+           
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Symmetric shape vs logCS - p-value = 0.004**")+ #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1) 


##New PCA plot with data corrected for allometry and symmetry
PCA_residuals_sym <- gm.prcomp(allometry_residuals_sym) 

#List of PC components and proportion of variations
PCA_residuals_sym

#Save PC scores as object to use later
pcscores_res_sym <- PCA_residuals_sym$x

#Save shapes of extremes for axes used in plot
PC1min_res_sym <- PCA_residuals_sym[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_res_sym <- PCA_residuals_sym[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_res_sym <- PCA_residuals_sym[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_res_sym <- PCA_residuals_sym[["shapes"]][["shapes.comp2"]][["max"]] 

#Show deformation grids on axis from mean shape, do this for all 4 extremes - "TPS" method
plotRefToTarget(mean_shape_residuals_sym, PC1min_res_sym, method = "TPS", mag = 1, label = FALSE)  #save image

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(mean_shape_residuals_sym, PC1min_res_sym, mesh = ref_mesh, method = "surface", mag = 1, label = FALSE)   #save as HTML

##3D windows save
#Save 3D window as html file - 3D widget
scene <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_res_sym <- as_tibble(pcscores_res_sym)
glimpse(pcscores_res_sym)
#Add labels and other attributes to tibble as columns
pcscores_res_sym <- pcscores_res_sym %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(pcscores_res_sym)

#Order tibble by variable e.g. age for plot legend
#Make factor for variable
pcscores_res_sym$age <- factor(pcscores_res_sym$age, 
                                levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
pcscores_res_sym <- pcscores_res_sym[order(pcscores_res_sym$age),]
glimpse(pcscores_res_sym)

#Nice plot
ggplot(pcscores_res_sym, aes(x = Comp1, y = Comp2, label = individuals, colour = age))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), 
                        values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+         #legend and color adjustments
  theme_bw()+
  xlab("PC 1 (47.38%)")+ #copy this from standard PCA plot or from PCA summary
  ylab("PC 2 (23.76%)")+
  ggtitle("PCA symmetry residuals")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position


#MODULARITY TEST ----

#Set all landmarks in one module
modules <- rep('a',16) 

#Put selected landmarks in second module
modules[4:9]<-'b' 
modules[14:16]<-'b' 
modules

#Perform modularity test - compare selected modules to the null assumption of random assignment of partitions (no modularity at all) 
#Is the data modular according to my defined modules?
#Best done on allometric residuals - NOT ON SYMMETRY data, not a biological significant hypothesis and little difference
modularity_test <- modularity.test(allometry_residuals, modules, CI = T,iter = 999, print.progress = T) 

#Get P value for CR values (same as RV)
summary(modularity_test) 

#Histogram of CR values
plot(modularity_test)

##Make better histogram plot with ggplot
#Save CR of data as object
CRdata <- modularity_test[["CR"]]

#Save all generated CRs as object
CRs <- modularity_test[["CR.boot"]]

#Calculate mean of CRs for arrow
mean_CR <- mean(CRs)

#Create tibble with modularity analysis values
modularity_plot_tibble <- data.frame(CRs)
modularity_plot_tibble <- as_tibble(modularity_plot_tibble)
glimpse(modularity_plot_tibble)

#Make simple histogram plot as object to obtain counts per bin with given bin width 
modularity_plot <- ggplot(modularity_plot_tibble, aes(CRs))+
  geom_histogram(binwidth = 0.01, fill = "gray91", colour = "gray78") #see if binwidth value appropriate and choose colors
modularity_plot

#Create data frame with plot variables
modularity_plot_data <- ggplot_build(modularity_plot)  

#Save only data counts as tibble
modularity_plot_data <- modularity_plot_data[["data"]][[1]]
modularity_plot_data <- modularity_plot_data[,1:5]
modularity_plot_data <- as_tibble(modularity_plot_data)

#Filter rows to select row with count for CR data and mean CR 
CRdata_filter <- modularity_plot_data %>% filter(xmin <= CRdata, xmax >= CRdata)
mean_CR_filter <- modularity_plot_data %>% filter(xmin <= mean_CR, xmax >= mean_CR)

#Create tibble with x and y columns to build arrows - x = mean position on bin, y = starting position on bin
modularity_plot_arrow <- data.frame(x_data = CRdata_filter$x, y_data = CRdata_filter$y, 
                                    x_mean = mean_CR_filter$x, y_mean = mean_CR_filter$y)
modularity_plot_arrow <- as_tibble(modularity_plot_arrow)

#Check that values of x are similar to the original data and mean values, if not change bins number or binwidth in original plot (add bins)
glimpse(modularity_plot_arrow)  

#Nice plot  
modularity_plot + #use plot obtained before after color and binwidth ok
 #add arrow for CR data
  geom_segment(data = modularity_plot_arrow, aes(x = x_data, xend = x_data, y = y_data, yend = y_data + 10, colour = "slateblue4"), size = 1.1,
               arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "first", type = "closed"), linejoin = 'mitre')+
 #add arrow for mean CR 
  geom_segment(data = modularity_plot_arrow, aes(x = x_mean, xend = x_mean, y = y_mean, yend = y_mean + 10, colour = "skyblue3"), size = 1.1,
               arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "first", type = "closed"), linejoin = 'mitre')+
#legend and color adjustments 
  scale_colour_manual(name = NULL, labels = c("mean CR (0.944)","CR data (0.882)"), #copy data from objects
                      values = c("slateblue4","skyblue3"))+            
  theme_minimal()+
  xlab("CR coefficient")+
  ylab("Frequency")+
  ggtitle("Modularity analysis - p-value = 0.003**")+  #copy data from standard plot
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))

#INTEGRATION TEST ----
##Perform integration test between the 2 modules (two-block PLS within configuration) 
#How much are the two or more modules integrated with each other?
#Best done on allometric residuals - NOT ON SYMMETRY data, not a biological significant hypothesis and little difference
integration_test <- integration.test(allometry_residuals, partition.gp = modules, iter = 999, print.progress = T) 

#Get P value and effect size Z - if significant means the modules are integrated!
summary(integration_test) 

#PLS plot of the two modules - two-block PLS with regression line
integration_pls_plot <- plot(integration_test, 
                        pch = 21, col = "darkgoldenrod2", bg = "darkgoldenrod2", cex = 1.2) 

#Save plot arguments as objects to use in plots
block1_pls_integr <- integration_pls_plot$plot.args$x
block2_pls_integr <- integration_pls_plot$plot.args$y

#Add labels
text(x = block1_pls_integr, y = block2_pls_integr, labels = specimens,
     pos = 3, offset = 0.5, cex = 0.75)   #improve appearance of labels


##Obtain shapes of modules at min and max of Block 1 and Block 2 axes

#Quick way to visualize shape changes on axes with interactive plot with "points" method
picknplot.shape(integration_pls_plot, label = TRUE, mag = 1)

#Better way using shape.predictor and plotReftoTarget - "TPS" and "surface" methods
#Save for block 1 and block 2 as objects from plot
block1_pls_integr_shapes <- integration_pls_plot$A1
block2_pls_integr_shapes <- integration_pls_plot$A2

#Create mean shape for each block
block1_pls_integr_mshape <- mshape(integration_pls_plot$A1)
block2_pls_integr_mshape <- mshape(integration_pls_plot$A2)

#Find specimen closer to mean for each block, useful to create warp mesh
findMeanSpec(block1_pls_integr_shapes)
findMeanSpec(block2_pls_integr_shapes)

#Create object containing only that specimen coordinates
block1_pls_integr_warp_spec <- block1_pls_integr_shapes[,,6]
block2_pls_integr_warp_spec <- block2_pls_integr_shapes[,,12]

#Import reference meshes for each block to use in surface method - parts of original reference mesh
block1_mesh_3D <- read.ply("Data/simpleskull_blocKs_tree1_rostrum.ply") #make sure NO binary encoding (ASCII)
block2_mesh_3D <- read.ply("Data/simpleskull_block2_braincase.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates to make sure it has same scale
range(block1_mesh_3D$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(block1_pls_integr_warp_spec)
range(block2_mesh_3D$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(block2_pls_integr_warp_spec)

#Create warp meshes, to use as reference for visualization of analyses
block1_pls_integr_ref_mesh <- warpRefMesh(mesh = block1_mesh_3D, mesh.coord = block1_pls_integr_warp_spec, 
                                          ref = block1_pls_integr_mshape, color = NULL, centered = FALSE) 
block2_pls_integr_ref_mesh <- warpRefMesh(mesh = block2_mesh_3D, mesh.coord = block2_pls_integr_warp_spec, 
                                          ref = block2_pls_integr_mshape, color = NULL, centered = FALSE) 

##Use shape predictor to generate shape coordinate at min and max of each block 
block1_pls_integr_shape_pred <- shape.predictor(block1_pls_integr_shapes, #shapes of block
                             x = block1_pls_integr, #values of axis form plot
                             Intercept = FALSE,  
                             method = "PLS",  #since it's the results of a PLS method
                             min = min(block1_pls_integr), #min value of axis
                             max = max(block1_pls_integr)) #max value of axis

#Show deformation grids on axis from mean shape, do this for all 4 extremes - "TPS" method
plotRefToTarget(block1_pls_integr_mshape, block1_pls_integr_shape_pred$min, mag = 1, method = "TPS") #save image
plotRefToTarget(block1_pls_integr_mshape, block1_pls_integr_shape_pred$max, mag = 1, method = "TPS") #save image

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(block1_pls_integr_mshape, block1_pls_integr_shape_pred$min, 
                mesh = block1_pls_integr_ref_mesh, method = "surface", mag = 1)   #save as HTML
plotRefToTarget(block1_pls_integr_mshape, block1_pls_integr_shape_pred$max, 
                mesh = block1_pls_integr_ref_mesh, method = "surface", mag = 1)   #save as HTML

##3D windows save
#Save screenshot of 3D window, useful for lateral and dorsal views - use screen snip if it fails
rgl.snapshot(filename = "Output/X.png") 
#Save 3D window as html file - 3D widget
blocKs_tree1_pred <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

#Repeat for other block
block2_pls_integr_shape_pred <- shape.predictor(block2_pls_integr_shapes, #shapes of block
                                     x = block2_pls_integr, #values of axis form plot
                                     Intercept = FALSE,  
                                     method = "PLS",  #since it's the results of a PLS method
                                     min = min(block2_pls_integr), #min value of axis
                                     max = max(block2_pls_integr)) #max value of axis

#Show deformation grids on axis from mean shape, do this for all 4 extremes - "TPS" method
plotRefToTarget(block2_pls_integr_mshape, block2_pls_integr_shape_pred$min, mag = 1, method = "TPS") #save image
plotRefToTarget(block2_pls_integr_mshape, block2_pls_integr_shape_pred$max, mag = 1, method = "TPS") #save image

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(block2_pls_integr_mshape, block2_pls_integr_shape_pred$min, 
                mesh = block2_pls_integr_ref_mesh, method = "surface", mag = 1)   #save as HTML
plotRefToTarget(block2_pls_integr_mshape, block2_pls_integr_shape_pred$max, 
                mesh = block2_pls_integr_ref_mesh, method = "surface", mag = 1)   #save as HTML

##3D windows save
#Save screenshot of 3D window, useful for lateral and dorsal views - use screen snip if it fails
rgl.snapshot(filename = "Output/X.png") 
#Save 3D window as html file - 3D widget
block2_pred <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)


##Make better integration PLS plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
integration_pls_plot_tibble <- data.frame(block1_pls_integr, block2_pls_integr)
integration_pls_plot_tibble

#Convert data frame to tibble
integration_pls_plot_tibble <- as_tibble(integration_pls_plot_tibble)
glimpse(integration_pls_plot_tibble)
#Add labels and other attributes to tibble as columns
integration_pls_plot_tibble <- integration_pls_plot_tibble %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(integration_pls_plot_tibble)

#Order tibble by variable e.g. age for plot legend
#Make factor for variable
integration_pls_plot_tibble$age <- factor(integration_pls_plot_tibble$age, 
                                levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
integration_pls_plot_tibble <- integration_pls_plot_tibble[order(integration_pls_plot_tibble$age),]
glimpse(integration_pls_plot_tibble)

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(integration_pls_plot_tibble, aes(x = block1_pls_integr, y = block2_pls_integr, label = individuals, colour = age))+
  #confidence intervals and reg line, before points  
  geom_smooth(method='lm',     #use standard function, values too small    
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+  
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), 
                        values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+          
  theme_classic(base_size = 12)+
  xlab("PLS1 Block 1: rostrum")+
  ylab("PLS1 Block 2: braincase")+
  ggtitle ("Integration test skull modules - p-value = 0.001***")+  #copy from test summary
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5), axis.title = element_text(size = 11))


#PHYLOMORPHOSPACE PLOT ----

#Import trees in Nexus format - branch lengths needed!!
trees1 <- "Data/trees_age_2.nex"   #trees with groups
trees2 <- "Data/trees_age.nex"     #specimen-level trees

##Read the trees for analysis
trees_age <- read.nexus(trees1) #trees with groups
plot(trees_age)
#Press enter to view next tree plots
trees_specimens <- read.nexus(trees2)   #specimen-level trees
plot(trees_specimens)
#Press enter to view next tree plots

#Check names of each tree in object
summary(trees_age)
summary(trees_specimens)    

##PCA Phylomorphospace
#Simple phylogeny projected on morphospace, no calculations
PCA_res_phylomor_tree1 <- gm.prcomp(allometry_residuals, phy = trees_specimens$age2)
PCA_res_phylomor_tree2 <- gm.prcomp(allometry_residuals, phy = trees_specimens$age4)

#View results
summary(PCA_res_phylomor_tree1)
summary(PCA_res_phylomor_tree2)

#Plot PCA with phylomorphospace
plot(PCA_res_phylomor_tree1, phylo = TRUE, #add phylogeny
     main = "PCA phylomorphospace - tree1",  pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 
plot(PCA_res_phylomor_tree2, phylo = TRUE, #add phylogeny
     main = "PCA phylomorphospace - tree2",  pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 

#3D phylo plot
plot(PCA_res_phylomor_tree1, time.plot = TRUE, pch = 21, col = "deeppink", bg = "deeppink", cex = 2, 
     phylo.par = list(edge.color = "grey60", edge.width = 1.5, tip.txt.cex = 0.75,
                      node.labels = F, anc.states = F))
plot(PCA_res_phylomor_tree2, time.plot = TRUE, pch = 21, col = "deeppink", bg = "deeppink", cex = 2, 
     phylo.par = list(edge.color = "grey60", edge.width = 1.5, tip.txt.cex = 0.75,
                      node.labels = F, anc.states = F))

##3D windows save
#Save screenshot of 3D window, useful for lateral and dorsal views - use screen snip if it fails
rgl.snapshot(filename = "Output/X.png") 
#Save 3D window as html file - 3D widget
PCA_3D <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

#PHYGENETIC SIGNAL TEST AND PaCa PLOT ----

##Test for phylogenetic signal in shape residuals
phylo_signal_res_tree1 <- physignal(allometry_residuals,phy = trees_specimens$age2, iter = 999)
phylo_signal_res_tree2 <- physignal(allometry_residuals,phy = trees_specimens$age4, iter = 999)

#View results
summary(phylo_signal_res_tree1)
summary(phylo_signal_res_tree2)

#Histogram of phylogenetic signal
plot(phylo_signal_res_tree1)
plot(phylo_signal_res_tree2)

##Make better histogram plot with ggplot (phylo_signal_res_tree1 as example)
#Save phylo signal K of data as object
Kdata_tree1 <- phylo_signal_res_tree1[["phy.signal"]]

#Save random Ks phylo signals as object
Ks_tree1 <- phylo_signal_res_tree1[["random.K"]]

#Calculate mean of Ks for arrow
mean_K_tree1 <- mean(Ks_tree1)

#Create tibble with phylogentic signal analysis values
phylo_signal_res_tree1_plot_tibble <- data.frame(Ks_tree1)
phylo_signal_res_tree1_plot_tibble <- as_tibble(phylo_signal_res_tree1_plot_tibble)
glimpse(phylo_signal_res_tree1_plot_tibble)

#Make simple histogram plot as object to obtain counts per bin with given bin width
phylo_signal_res_tree1_plot <- ggplot(phylo_signal_res_tree1_plot_tibble, aes(Ks_tree1))+
  geom_histogram(bins = 30, fill = "gray91", colour = "gray78")     #see if bins value appropriate and choose colors
phylo_signal_res_tree1_plot

#Create data frame with plot variables
phylo_signal_res_tree1_plot_data <- ggplot_build(phylo_signal_res_tree1_plot)  

#Save only data counts as tibble
phylo_signal_res_tree1_plot_data <- phylo_signal_res_tree1_plot_data[["data"]][[1]]
phylo_signal_res_tree1_plot_data <- phylo_signal_res_tree1_plot_data[,1:5]
phylo_signal_res_tree1_plot_data <- as_tibble(phylo_signal_res_tree1_plot_data)
glimpse(phylo_signal_res_tree1_plot_data)

#Filter rows to select row with count for CR data and mean CR 
Kdata_tree1_filter <- phylo_signal_res_tree1_plot_data %>% filter(xmin <= Kdata_tree1, xmax >= Kdata_tree1)
mean_K_tree1_filter <- phylo_signal_res_tree1_plot_data %>% filter(xmin <= mean_K_tree1, xmax >= mean_K_tree1)

#Create tibble with x and y columns to build arrows - x = mean position on bin, y = starting position on bin
phylo_signal_res_tree1_arrow_plot <- data.frame(x_data = Kdata_tree1_filter$x, y_data = Kdata_tree1_filter$y, 
                                                x_mean = mean_K_tree1_filter$x, y_mean = mean_K_tree1_filter$y)
phylo_signal_res_tree1_arrow_plot <- as_tibble(phylo_signal_res_tree1_arrow_plot)

#Check that values of x are similar to the original data and mean values, if not change bins number or binwidth in original plot (add bins)
glimpse(phylo_signal_res_tree1_arrow_plot)

#Nice plot  
phylo_signal_res_tree1_plot + #use plot obtained before after color and binwidth ok
  #add arrow for CR data
  geom_segment(data = phylo_signal_res_tree1_arrow_plot, aes(x = x_data, xend = x_data, y = y_data, yend = y_data + 20, colour = "firebrick4"), size = 1,
               arrow = arrow(angle = 30, length = unit(0.02, "npc"), ends = "first", type = "closed"), linejoin = 'mitre')+
  #add arrow for mean CR 
  geom_segment(data = phylo_signal_res_tree1_arrow_plot, aes(x = x_mean, xend = x_mean, y = y_mean, yend = y_mean + 20, colour = "black"), size = 1,
               arrow = arrow(angle = 30, length = unit(0.02, "npc"), ends = "first", type = "closed"), linejoin = 'mitre')+
#legend and color adjustments  
  scale_colour_manual(name = NULL, labels = c("mean K (0.348)", "Observed K (0.325)"),    #copy data from objects
                      values = c("black", "firebrick4"))+        
  theme_minimal()+
  xlab("Phylogenetic Signal K")+
  ylab("Frequency")+
  ggtitle("Phylogentic Signal test tree1 - p-value = 0.743")+  #copy data from standard plot
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))

#Phylogenetic signal test for shapes produces PaCa
##PaCA: Phylogenetically-aligned PCA
#Projection of phylogenetic signal in the first two components, helps to highlight importance of phylogeny relative to other signals
PCA_res_PaCA_tree1 <- phylo_signal_res_tree1$PACA
PCA_res_PaCA_tree2 <- phylo_signal_res_tree2$PACA

#View results
summary(PCA_res_PaCA_tree1)
summary(PCA_res_PaCA_tree2)

#Plot PaCa
plot(PCA_res_PaCA_tree1, phylo = TRUE, #add phylogeny
     main = "PaCa tree1", pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 
plot(PCA_res_PaCA_tree2, phylo = TRUE, #add phylogeny
     main = "PaCa tree2", pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2)

#phyloPCA PLOT ----
#Phylogeny is used to calculate principal components of variation - GLS method

#PCA analysis with phylogeny
PCA_res_phylo_tree1 <- gm.prcomp(allometry_residuals, phy = trees_specimens$age2, GLS = TRUE)
PCA_res_phylo_tree2 <- gm.prcomp(allometry_residuals, phy = trees_specimens$age4, GLS = TRUE)

#View results
summary(PCA_res_phylo_tree1)
summary(PCA_res_phylo_tree2)

#Plot PCA
plot(PCA_res_phylo_tree1, phylo = TRUE, #add phylogeny
     main = "phyloPCA tree1", pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 
plot(PCA_res_phylo_tree2, phylo = TRUE, #add phylogeny
     main = "phyloPCA tree2", pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 

##Make better PCA plot using ggplot - phyloPCA
#Tree 1 as example
#Save pc scores as object
pcscores_res_phylo_tree1 <- PCA_res_phylo_tree1[["x"]]

#Read PC scores as tibble
pcscores_res_phylo_tree1  <- as_tibble(pcscores_res_phylo_tree1)
glimpse(pcscores_res_phylo_tree1)
#Add labels and other attributes to tibble as columns
pcscores_res_phylo_tree1  <- pcscores_res_phylo_tree1  %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(pcscores_res_phylo_tree1)

##Order tibble by variable e.g. age
#Make factor for variable
pcscores_res_phylo_tree1$age <- factor(pcscores_res_phylo_tree1$age, 
                                       levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
pcscores_res_phylo_tree1 <- pcscores_res_phylo_tree1[order(pcscores_res_phylo_tree1$age),]
glimpse(pcscores_res_phylo_tree1)

#Use package ggphylomorpho to create a first phylomorphospace plot
PCA_res_phylo_tree1_plot <- ggphylomorpho(tree = trees_specimens$age2, tipinfo = pcscores_res_phylo_tree1, 
                                xvar = Comp1, yvar = Comp2, 
                                factorvar = age, labelvar = individuals, tree.alpha = 0.7)
PCA_res_phylo_tree1_plot 

#Using the package gginnards, remove the GeomPoint (points on tips) and GeomTextRepel (labels) so you just have the tree network
PCA_res_phylo_tree1_plot  <- PCA_res_phylo_tree1_plot  %>%
  delete_layers("GeomPoint") %>%   delete_layers("GeomTextRepel")

#Create convex hulls - polygons that include all specimens that belong to a group
hulls_tree1 <- pcscores_res_phylo_tree1 %>%
  group_by(age) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)
glimpse(hulls_tree1)

#Nice plot with phylogeny + convex hulls
PCA_res_phylo_tree1_plot  + 
  geom_point(data = pcscores_res_phylo_tree1, aes(x = Comp1, y = Comp2, colour = age), size = 3)+ #add "shape =" to aes to change shapes based on groups
#Add convex hulls
  geom_polygon(data = hulls_tree1, aes(x = x, y = y, fill = age), alpha = .5, show.legend = FALSE)+ #alpha = transparency of colour
#Add labels  
  geom_text_repel(data = pcscores_res_phylo_tree1, aes(x = Comp1, y = Comp2, label = individuals, size = 3.5), show.legend = FALSE)+
#Colour points and convex hulls
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+            #legend and color adjustments
  scale_fill_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"),
                    values = c("cyan2","deepskyblue1","dodgerblue3","blue4"))+
  theme_bw()+
  ggtitle("phyloPCA tree1")+
  theme(legend.title = element_text(face="bold"), #Legend titles in bold
        legend.justification = "top")+ #Legend position
  xlab("PC 1 (41.92%)")+ #copy this from standard PCA plot
  ylab("PC 2 (27.47%)")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position

##Extra code ggplot
#Code to remove dots from legend
  guides(fill = guide_legend(override.aes = list(shape = NA)))+  #remove annoying dots in the colour legend
#Code to decide which polygons are sued for each group if shape = present
scale_shape_manual(values=c(17, 18, 15)) 
  

#ALLOMETRY ANALYSIS BY GROUP AND GROUP MEANS ----

#Regression of shapes on logCS considering groups as additional factor (e.g. species, age) - test if grouping is driving allometry
allometry_group <- procD.lm(coords ~ logCsize * factor_age, iter=999, print.progress = TRUE) 
View(allometry_group)

#Main results of ANOVA analysis of allometry with logCS
summary(allometry_group) 

##Plot shape vs logCS to visualize allometry
#Diagnostic plots to check if model is appropriate - similar to ANOVA tables
init <- par(no.readonly=TRUE) #store initial plot parameters to restore later
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
plot(allometry_group,type = "diagnostics", cex = 1.2, font.main = 2)
par(init)                     #restore initial plot parameters (1 plot showing at a time)

#Regression score of shape vs logCS - regression method with regression score plotting
allometry_group_plot_regscore <- plot(allometry_group, type = "regression", predictor = logCsize, reg.type = "RegScore",
                                main = "Shape vs logCS by group",xlab = "logCS", 
                                pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics
text(x = logCsize, y = allometry_group_plot_regscore$RegScore, labels = specimens,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

##Add regression line with confidence intervals to plot
#Create object to use for linear model
allometry_group_regscores <- allometry_group_plot_regscore[["RegScore"]] 

#Linear model for line
allometry_group_regline <- (lm(allometry_group_regscores ~ logCsize))

#Check p-value for plot
summary(allometry_group_regline)

#Draw line on plot
abline(allometry_group_regline, col = "darkseagreen3", 
       lty = 2, lwd = 2)     #line type (e.g. dashed) and width

#Add confidence intervals
#Create data for confidence intervals
allometry_group_xvals <- seq(min(logCsize), max(logCsize), length = 12)   #use min and max of x values (logCS) as limits and use number of specimens as length of sequence
allometry_group_newX <- expand.grid(logCsize = allometry_group_xvals)     #warp x_vals on values of x axis (logCS)
allometry_group_newY <- predict(allometry_group_regline, newdata = data.frame(x = allometry_group_newX), 
                                interval="confidence", level = 0.95)      #predict the y values based on the x sequence

#Draw confidence intervals lines on plot
matlines(allometry_group_newX, allometry_group_newY[,2:3],  #first column of newY not useful, it is the fit, 2 and 3 are the min and max values
         col = "darkseagreen3", lty=1)                     #line graphics

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_group_plot_tibble <- data.frame(logCS = allometry_group_plot_regscore[["plot.args"]][["x"]], 
                                          RegScores = allometry_group_plot_regscore[["plot.args"]][["y"]])
allometry_group_plot_tibble

#Convert data frame to tibble
allometry_group_plot_tibble <- as_tibble(allometry_group_plot_tibble)
glimpse(allometry_group_plot_tibble)
#Add labels and other attributes to tibble as columns
allometry_group_plot_tibble <- allometry_group_plot_tibble %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(allometry_group_plot_tibble)

#Order tibble by variable e.g. age for plot legend
#Make factor for variable
allometry_group_plot_tibble$age <- factor(allometry_group_plot_tibble$age, 
                                levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
allometry_group_plot_tibble <- allometry_group_plot_tibble[order(allometry_group_plot_tibble$age),]
glimpse(allometry_group_plot_tibble)

##Add regression line with confidence intervals
#Make data frame of data for confidence intervals
allometry_group_conf_intervals <- data.frame(allometry_group_newX, allometry_group_newY)
allometry_group_conf_intervals 
#Rename columns to match main plot tibble varibales for x and y
allometry_group_conf_intervals <- rename(allometry_group_conf_intervals, logCS = logCsize, RegScores = fit)
allometry_group_conf_intervals 

#Convert data frame to tibble
allometry_group_conf_intervals <- as_tibble(allometry_group_conf_intervals)
glimpse(allometry_group_conf_intervals)
#Add labels and other attributes to tibble as columns to match main plot tibble
allometry_group_conf_intervals <- allometry_group_conf_intervals %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
glimpse(allometry_group_conf_intervals)

#Order tibble by variable e.g. age for plot legend
#Make factor for variable
allometry_group_conf_intervals$age <- factor(allometry_group_conf_intervals$age, 
                                          levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
allometry_group_conf_intervals <- allometry_group_conf_intervals[order(allometry_group_conf_intervals$age),]
glimpse(allometry_group_conf_intervals)

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(allometry_group_plot_tibble, aes(x = logCS, y = RegScores, label = individuals, colour = age))+
  geom_smooth(data = allometry_group_conf_intervals, aes(ymin = lwr, ymax = upr), stat = 'identity',  #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), 
                        values = c("cyan2","deepskyblue1","dodgerblue3", "blue4"))+          
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Shape vs logCS by group - p-value = 0.765")+  #copy p-value from linear model of conf_intervals, gives idea if groups are significant factor or not
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1)     #position of tables relative to point (proximity and distance) 

##Group means - can be useful if I have multiple specimens of 1 group to plot as a single point
#Create mean shape for each group based on classifiers
shape_group_means <- aggregate(two.d.array(coords) ~ factor_age, FUN = mean)
shape_group_means

#Set row names that match tree labels
rownames(shape_group_means) <- c("adult","earlyFetus","lateFetus","neonate") #check that the factors are in the same order as the rows currently
shape_group_means

#Transform in matrix and exclude columns that do not contain coordinates
shape_group_means <- as.matrix(shape_group_means[,-1]) #eliminate first column numerically
#Transform in 3D array
shape_group_means <- arrayspecs(shape_group_means,16,3) 

#Mean logCS size based on groups
size_group_means <- aggregate(logCsize ~ factor_age, FUN = mean) 

#Set row names that match tree labels
rownames(size_group_means) <- c("adult","earlyFetus","lateFetus","neonate") 

#Eliminate first column with names of classifiers
size_group_means <- select(size_group_means,-"factor_age")#eliminate first column by selecting the rest and writing column name

#Make numeric
size_group_means <- as.matrix(size_group_means) 
