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


#DATA IMPORT AND PREP ----

#Import data into R, make sure only first columns are not numbers
minke_raw <- read_csv("Data/16minkeonly.csv") 
glimpse(minke_raw)

#Average landmark takes
minke_avg_takes <- minke_raw %>% group_by(specimenID) %>% summarize(across(starts_with("Raw"), list(mean)))
glimpse(minke_avg_takes)

#Create data frame with only numerical values and set row names as specimen ID
minke_avg_takes <- remove_rownames(minke_avg_takes)
has_rownames(minke_avg_takes) #steps for tibbles
minke_avg_takes2 <- data.frame(column_to_rownames(minke_avg_takes, var = "specimenID")) #adding row names to data frame
minke_avg_takes <- minke_avg_takes2 
remove(minke_avg_takes2) #change name and remove extra object

#Save individual names as object
minke_ind <- row.names(minke_avg_takes) 

#Transform in 3D array, first number is number of landmarks, second is dimensions (2 or 3)
minke_array <- arrayspecs(minke_avg_takes, 16, 3) 

##Extract classifier columns from raw data
#Make tibble for plots in ggplot
minke_age_tibble <- minke_raw %>% 
  group_by(specimenID) %>% 
  summarize(age) %>%       #this creates a tibble with only the individual names and one text factor (e.g. age)
  distinct()               #this only keeps in the tibble unique combinations of the two values, so you get a tibble with the same number or rows as the data frame
glimpse(minke_age_tibble)
has_rownames(minke_age_tibble)    #always check

#Save classifiers as factors from tibble
minke_age <- as.factor(minke_age_tibble$age)

#remove objects that are not needed (e.g. "links")
remove(A) 


#GPA ALIGNMENT ----

#Procrustes alignment, should also show mean config coordinates
minke_gpa<-gpagen(minke_array) 
View(minke_gpa)
plot(minke_gpa) #see points in space

#Save Centroid size as object
Csize<-minke_gpa$Csize 
#Log-transform Centroid size as object
logCsize<- log10(Csize) 

#Save mean shape to create links
mean_shape<-minke_gpa$consensus 

#Coordinates of all specimens after GPA alignment
minke_coords<-minke_gpa$coords 

#PREPARE WARP MESH AND LINKS  ----

#Find specimen closer to mean, useful to create warp mesh
findMeanSpec(minke_coords) #number below specimen name is the number of the specimen in the array

#Create object containing only that specimen coordinates
mean_spec <- minke_coords[,,2]
mean_spec

#Import simplified mesh to create warp mesh on
minke_skull3D <- read.ply("Data/simpleskull.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates to make sure it has same scale
range(minke_skull3D$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(mean_spec)

##Create warp mesh, to use as reference for visualization of analyses
ref_mesh <- warpRefMesh(mesh = minke_skull3D, mesh.coord = mean_spec, ref = mean_shape, color = NULL, centered = FALSE) 

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
PCA_all<-gm.prcomp(minke_coords)

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
minke_pcscores_all <- PCA_all$x 

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
PC1min <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

##Make better PCA plot using ggplot
#Read PC scores as tibble
minke_pcscores_all <- as_tibble(minke_pcscores_all)
glimpse(minke_pcscores_all)
View(minke_pcscores_all)
#Add labels and other attributes to tibble as columns
minke_pcscores_all <- minke_pcscores_all %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(minke_pcscores_all)

#Nice plot
ggplot(minke_pcscores_all, aes(x = Comp1, y = Comp2, label = individuals, colour = age))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+            #legend and color adjustments
  theme_bw()+
  xlab("PC 1 (48.89%)")+ #copy this from standard PCA plot
  ylab("PC 2 (22.17%)")+
  ggtitle("PCA all data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position

##Regression PC1 and PC2 vs logCS
#Create data frame with data
minke_pcscores_all_df <- data.frame(PCA_all$x , size = logCsize)
minke_pcscores_all_df

#Calculate regression for each component
PC1_size_reg <- lm(Comp1~size, data = minke_pcscores_all_df)
PC2_size_reg <- lm(Comp2~size, data = minke_pcscores_all_df)

#View results and p-value
summary(PC1_size_reg)
summary(PC2_size_reg)

#Diagnostic plots for both regressions
init <- par(no.readonly=TRUE) #store initial plot parameters to restore later
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
plot(PC1_size_reg, cex = 1.2, font.main = 2)
plot(PC2_size_reg, cex = 1.2, font.main = 2)                             
par(init)                     #restore initial plot parameters (1 plot showing at a time)

#Plot regression with ggplot
#Convert data frame to tibble
minke_pcscores_all_tibble <- as_tibble(minke_pcscores_all_df)
glimpse(minke_pcscores_all_tibble)
#Add labels and other attributes to tibble as columns
minke_pcscores_all_tibble <- minke_pcscores_all_tibble %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(minke_pcscores_all_tibble)

#Create data frame with line parameters from regression
#Allows to show a line on PC plot with specimens colored IF col and other graphics OUTSIDE of aes()!!!
PC_reg_line <- data.frame(int1 = PC1_size_reg[["coefficients"]][["(Intercept)"]], slope1 = PC1_size_reg[["coefficients"]][["size"]], #PC1 values
                          int2 = PC2_size_reg[["coefficients"]][["(Intercept)"]], slope2 = PC2_size_reg[["coefficients"]][["size"]]) #PC2 values
PC_reg_line

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(minke_pcscores_all_tibble, aes(x = size, y = Comp1, label = individuals, colour = age))+
  #line on plot
  geom_abline(data = PC_reg_line, aes(intercept = int1, slope = slope1), colour = "darkblue", size = 0.8, linetype = "dashed", show.legend = FALSE)+
  geom_point(size = 3)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  xlab("logCS")+
  ylab("PC 1 (48.89%)")+
  geom_text_repel(colour = "black", size = 3.5)

#Repeat for other component
#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(minke_pcscores_all_tibble, aes(x = size, y = Comp2, label = individuals, colour = age))+
  #line on plot
  geom_abline(data = PC_reg_line, aes(intercept = int2, slope = slope2), colour = "darkblue", size = 0.8, linetype = "dashed", show.legend = FALSE)+
  geom_point(size = 3)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  xlab("logCS")+
  ylab("PC 2 (22.17%)")+
  geom_text_repel(colour = "black", size = 3.5)


#ALLOMETRY CORRECTION ----
##Evaluate allometry and get the allometry-free shapes using LogCS, use this for analyses
minkeAllometry_log <- procD.lm(minke_coords~logCsize, iter=999, print.progress = TRUE) 
View(minkeAllometry_log)

#Main results of ANOVA analysis of allometry with logCS
summary(minkeAllometry_log) 

#Create residuals array to then save as coordinates for analyses
shape_residuals <- arrayspecs(minkeAllometry_log$residuals,p=dim(minke_coords)[1], k=dim(minke_coords)[2]) 

#New shapes adjusted for allometry with CS to use in analyses
minkeAllometry_residuals <- shape_residuals + array(mean_shape, dim(shape_residuals)) 
View(minkeAllometry_residuals)

#Save mean shape of allometry-adjusted shapes to sue later
mean_shape_res <- mshape(minkeAllometry_residuals)

##Plot shape vs logCS to visualize allometry
#Diagnostic plots to check if model is appropriate - similar to ANOVA tables
init <- par(no.readonly=TRUE) #store initial plot parameters to restore later
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
allometryplot_diagnostics <- plot(minkeAllometry_log,type = "diagnostics",
                              cex = 1.2, font.main = 2)
par(init)                     #restore initial plot parameters (1 plot showing at a time)

#Regression score of shape vs logCS - regression method with regression score plotting
allometryplot_regscore <- plot(minkeAllometry_log,type = "regression",predictor = logCsize, reg.type = "RegScore",
                          main = "Shape vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics
text(x = logCsize, y = allometryplot_regscore$RegScore, labels = minke_ind,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

##Add regression line with confidence intervals to plot
#Create object to use for linear model
reg_scores <- allometryplot_regscore[["RegScore"]] 

#Linear model for line
reg_line <- (lm(reg_scores~logCsize))

#Draw line on plot
abline(reg_line, col = "darkseagreen3", 
       lty = 2, lwd = 2)     #line type (e.g. dashed) and width

#Add confidence intervals
#Create data for confidence intervals
x_vals <- seq(min(logCsize), max(logCsize), length = 12)   #use min and max of x values (logCS) as limits and use number of specimens as length of sequence
newX <- expand.grid(logCsize = x_vals)                     #warp x_vals on values of x axis (logCS)
newY <- predict(reg_line, newdata = data.frame(x = newX), interval="confidence",
        level = 0.95)                                      #predict the y values based on the x sequence

#Draw confidence intervals lines on plot
matlines(newX, newY[,2:3],                                 #first column of newY not useful, it is the fit, 2 and 3 are the min and max values
         col = "darkseagreen3", lty=1)                     #line graphics

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
minkeAllometry_plot <- data.frame(logCS = allometryplot_regscore[["plot.args"]][["x"]], RegScores = allometryplot_regscore[["plot.args"]][["y"]])
minkeAllometry_plot

#Convert data frame to tibble
minkeAllometry_plot <- as_tibble(minkeAllometry_plot)
glimpse(minkeAllometry_plot)
View(minkeAllometry_plot)
#Add labels and other attributes to tibble as columns
minkeAllometry_plot <- minkeAllometry_plot %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(minkeAllometry_plot)

##Add regression line with confidence intervals
#Make data frame of data for confidence intervals
conf_intervals <- data.frame(newX, newY)
conf_intervals 
#Rename columns to match main plot tibble varibales for x and y
conf_intervals <- rename(conf_intervals, logCS = logCsize, RegScores = fit)
conf_intervals 

#Convert data frame to tibble
conf_intervals <- as_tibble(conf_intervals)
glimpse(conf_intervals)
#Add labels and other attributes to tibble as columns to match main plot tibble
conf_intervals <- conf_intervals %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(conf_intervals)

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(minkeAllometry_plot, aes(x = logCS, y = RegScores, label = individuals, colour = age))+
  geom_smooth(data = conf_intervals, aes(ymin = lwr, ymax = upr), stat = 'identity',          #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Shape vs logCS")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1)     #position of tables relative to point (proximity and distance) 


#PC1 values vs logCS - regression method with prediction line plotting
allometryplot_regline <- plot(minkeAllometry_log,type = "regression",predictor = logCsize, reg.type = "PredLine",
                             main = "PC1 vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)
text(x = logCsize, y = allometryplot_regline$PredLine, labels = minke_ind,
     pos = 2, offset = 0.5, cex = 0.75)  

#PCA plot of fitted values - likely not very useful if size is a big component of variation, only PC1 will have weight
allometryplot_pc <- plot(minkeAllometry_log,type = "PC",
                         main = "PCA fitted values", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)
text(x = allometryplot_pc$plot.args$x, y = allometryplot_pc$plot.args$y, labels = minke_ind,
     pos = 3, offset = 0.5, cex = 0.75) 

##Specific plotAllometry functions - "CAC" and "size.shape" PCA (unclear what is the purpose of size-shape PCA, will not use)
#CAC plot - for simple allometry is the same as RegScore plot and two-block PLS plot of shape and size
allometryplot_CAC <- plotAllometry(minkeAllometry_log, size = logCsize, logsz = FALSE, method = "CAC", 
                      main = "CAC plot", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics

##Make better CAC allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
minkeAllometry_CAC_plot <- data.frame(logCS = allometryplot_CAC[["size.var"]], CAC = allometryplot_CAC[["CAC"]], #common allometric component
                                      RSC1 = allometryplot_CAC[["all.plot.args"]][["RSC"]][["x"]])               #main residual shape component
minkeAllometry_CAC_plot

#Convert data frame to tibble
minkeAllometry_CAC_plot <- as_tibble(minkeAllometry_CAC_plot)
glimpse(minkeAllometry_CAC_plot)
#Add labels and other attributes to tibble as columns
minkeAllometry_CAC_plot <- minkeAllometry_CAC_plot %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(minkeAllometry_CAC_plot)

##Two plots: CAC vs logCsize (A) and RSC1 vs CAC (B)
#Plot A
##Add regression line with confidence intervals
#Create object to use for linear model
reg_scores_CAC <- allometryplot_CAC[["CAC"]]

#Linear model for line
reg_line_CAC <- (lm(reg_scores_CAC~logCsize))

#Add confidence intervals
#Create data for confidence intervals
x_vals <- seq(min(logCsize), max(logCsize), length = 12)   #use min and max of x values (logCS) as limits and use number of specimens as length of sequence
newX <- expand.grid(logCsize = x_vals)                     #warp x_vals on values of x axis (logCS)
newY <- predict(reg_line_CAC, newdata = data.frame(x = newX), interval="confidence",
                level = 0.95)                                      #predict the y values based on the x sequence

#Make data frame of data for confidence intervals
conf_intervals_CAC <- data.frame(newX, newY)
#Rename columns to match main plot tibble varibales for x and y
conf_intervals_CAC <- rename(conf_intervals_CAC, logCS = logCsize, CAC = fit)
conf_intervals_CAC

#Convert data frame to tibble
conf_intervals_CAC <- as_tibble(conf_intervals_CAC)
glimpse(conf_intervals_CAC)
#Add labels and other attributes to tibble as columns to match main plot tibble
conf_intervals_CAC <- conf_intervals_CAC %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(conf_intervals_CAC)

#Nice plot with specimens colored by age AND regression line with confidence intervals  - no labels
ggplot(minkeAllometry_CAC_plot, aes(x = logCS, y = CAC, label = individuals, colour = age))+
  geom_smooth(data = conf_intervals_CAC, aes(ymin = lwr, ymax = upr), stat = 'identity',          #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  ggtitle ("CAC vs logCS")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Plot B

#Nice plot with specimens colored by age - no labels
ggplot(minkeAllometry_CAC_plot, aes(x = CAC, y = RSC1, colour = age))+
 #confidence intervals and reg line using standard function, difficult to do with external model for residuals (too small)
 #see function in CAC plot A for code for external model
  geom_smooth(method = 'lm',        #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.5)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  ggtitle ("Residual shape component (RSC) 1 vs CAC")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


#TWO-BLOCK PLS ----
#Two-block PLS of allometry, another way to visualize connection between logCS and shape, main one for analyses
minkeAllometry_pls <- two.b.pls(logCsize, minke_coords,iter = 999) 

#Get P-value of regression
minkeAllometry_pls

#Plot two-block PLS with regression line
allometryplot_pls <- plot(minkeAllometry_pls, 
      pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2)   #improve appearance of points
#Save plot arguments as objects to use in plots
block1 <- allometryplot_pls$plot.args$x
block2 <- allometryplot_pls$plot.args$y
#Add labels
text(x = block1, y = block2, labels = minke_ind,
     pos = 3, offset = 0.5, cex = 0.75)   #improve appearance of labels

##Make better PLS plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
minkeAllometry_pls_plot <- data.frame(block1, block2)
minkeAllometry_pls_plot

#Convert data frame to tibble
minkeAllometry_pls_plot <- as_tibble(minkeAllometry_pls_plot)
glimpse(minkeAllometry_pls_plot)
View(minkeAllometry_pls_plot)
#Add labels and other attributes to tibble as columns
minkeAllometry_pls_plot <- minkeAllometry_pls_plot %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(minkeAllometry_pls_plot)

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(minkeAllometry_pls_plot, aes(x = block1, y = block2, label = individuals, colour = age))+
#confidence intervals and reg line, before points  
  geom_smooth(method='lm',     #use standard function, values too small    
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+  
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  xlab("PLS1 Block 1: logCS")+
  ylab("PLS1 Block 1: Shape")+
  ggtitle ("PLS1 plot: Block 1 (logCS) vs Block 2 (Shape)")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5), axis.title = element_text(size = 11))


#PCA ALLOMETRY RESIDUALS ----
#New PCA plot with data corrected for allometry
PCA_residuals <- gm.prcomp(minkeAllometry_residuals) 

#List of PC components and proportion of variations
PCA_residuals

##View plot
plot(PCA_residuals, main = "PCA residuals",  pch = 21, #title and type of point to be used
     col = "deeppink", bg = "deeppink", cex = 1, font.main = 2)  #improve graphics
#Add quick labels to plot
text(x = PCA_residuals$x[,1], y = PCA_residuals$x[,2], labels = rownames(PCA_residuals$x), 
     pos = 1, offset = 0.5, cex = 0.75)    #improve graphics

#Save PC scores as object to use later
minke_pcscores_res <- PCA_residuals$x

#Save shapes of extremes for axes used in plot
PC1min_res <- PCA_residuals[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_res <- PCA_residuals[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_res <- PCA_residuals[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_res <- PCA_residuals[["shapes"]][["shapes.comp2"]][["max"]] 

#Show deformation grids on axis from mean shape, do this for all 4 extremes - "TPS" method
plotRefToTarget(mean_shape_res, PC1min_res, method = "TPS", mag = 1, label = FALSE)  #save image

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - "points" method
plotRefToTarget(mean_shape_res, PC1min_res, method = "points", mag = 1, label = FALSE)   #save as HTML

#Show 3D deformation from mean with vectors, do this for all 4 extremes - "vector" method
plotRefToTarget(mean_shape_res, PC1min_res, method = "vector", mag = 1, label = FALSE)   #save as screenshot

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(mean_shape_res, PC1min_res, mesh = ref_mesh, method = "surface", mag = 1, label = FALSE)   #save as HTML

##3D windows save
#Save screenshot of 3D window, useful for lateral and dorsal views - use screen snip if it fails
rgl.snapshot(filename = "Output/X.png") 
#Save 3D window as html file - 3D widget
PC1min <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

##Make better PCA plot using ggplot
#Read PC scores as tibble
minke_pcscores_res <- as_tibble(minke_pcscores_res)
glimpse(minke_pcscores_res)
View(minke_pcscores_res)
#Add labels and other attributes to tibble as columns
minke_pcscores_res <- minke_pcscores_res %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(minke_pcscores_res)

#Nice plot
ggplot(minke_pcscores_res, aes(x = Comp1, y = Comp2, label = individuals, colour = age))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+            #legend and color adjustments
  theme_bw()+
  xlab("PC 1 (44.15%)")+ #copy this from standard PCA plot
  ylab("PC 2 (22.79%)")+
  ggtitle("PCA residuals")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position


#ANOVA and PROCUSTES VARIANCES OF GROUP DIFFERENCES  ----
#Create dataframe to operate more easily - use allometry residuals
minke_dataframe <- geomorph.data.frame(minkeAllometry_residuals, gp = minke_age) 

#Conduct ANOVA to test differences between groups
minke_age_anova <- procD.lm(minkeAllometry_residuals ~ gp,iter=999, data = minke_dataframe) 

#Results and significance of ANOVA
summary(minke_age_anova) 

#Diagnostic plots to check if model is appropriate - similar to ANOVA tables - DO also if not significant, don't do other plots
init <- par(no.readonly=TRUE) #store initial plot parameters to restore later
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
age_anova_diagnostics <- plot(minke_age_anova,type = "diagnostics",
                                  cex = 1.2, font.main = 2)
par(init)                     #restore initial plot parameters (1 plot showing at a time)


##Morphological disparity 
#Calculate Procrustes variances and distances between groups after allometric correction
minke_age_disparity <- morphol.disparity(minkeAllometry_residuals ~ 1, groups = ~ minke_age, iter = 999) #do NOT use dataframe, issues with groups

#Results and significance
summary(minke_age_disparity) 

#TRAJECTORY ANALYSIS ----
#Shows trajectories of variation using groups, use obj from procD.lm
#Trajectory points must be defined and they can either be a factor to assess different trajectories within groups or they should be = to the nummber of groups
trajectory_age <- trajectory.analysis(minkeAllometry_log, groups = minke_age_tibble$age, traj.pts = 2, pca = TRUE, print.progress = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different?
summary(trajectory_age, show.trajectories = TRUE, attribute = "MD") 
#Trajectory correlations - only useful if there are significant differences 
summary(trajectory_age, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
#Trajectory shape differences - can only be used if each trajectory/group has 3 or more points/specimens - what is the distance between trajectories?
summary(trajectory_age, show.trajectories = TRUE, attribute = "SD") 

#Plot results - PCA of fitted values, similar to CVA
trajectory_plot <- plot(trajectory_age, main = "Trajectories by age",  pch = 21, #title and type of point to be used
                        col = "gray", bg = "gray", cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_plot, traj.pch = 21, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, start.bg = 3, end.bg = 2) 

#doesn't work -> legend("right", legend = minke_age)

##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_pcscores <- as_tibble(trajectory_plot[["pc.points"]])
glimpse(trajectory_pcscores)

#Find out order of group variables
glimpse(trajectory_plot[["trajectories"]])

#Create data frame that contains the group variables in order and that has = number of rows to pc score points
rows_to_repeat <- bind_rows(data.frame(group = "adult", counter = 1:12), data.frame(group = "earlyFetus", counter=13:24), 
                            data.frame(group ="lateFetus", counter=25:36), data.frame(group = "neonate", counter=37:48))
glimpse(rows_to_repeat)
#Delete extra counter column
rows_to_repeat <- rows_to_repeat [,-2]

#Add group names and other attributes to tibble as columns
trajectory_pcscores <- trajectory_pcscores %>% mutate(age = rows_to_repeat)
glimpse(trajectory_pcscores)

#Calculate means of PC1 and PC2 per group to draw trajectories
trajectory_pcscores_means <- trajectory_pcscores %>% group_by(age)  %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_pcscores_means <- trajectory_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
trajectory_pcscores_means

#Nice plot
ggplot(trajectory_pcscores, aes(x = PC1, y = PC2, colour = age))+
  geom_point(size = 3)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+            #legend and color adjustments
  theme_bw()+
  xlab("PC 1 (60.17%)")+ #copy this from standard PCA plot
  ylab("PC 2 (30.92%)")+
  ggtitle("Trajectories by age")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+  #title font and position
 #add trajectory lines, one line for each, write column name and row number form tibble
  geom_segment(data = trajectory_pcscores_means, aes(x = x[4], y = y[4], #neonate
                   xend = x[1], yend = y[1]),  #adult 
                   colour = "gray35", size = 0.8, linejoin = 'mitre', arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))+  #add arrow at end
  geom_segment(data = trajectory_pcscores_means, aes(x = x[3], y = y[3], #lateF
                                                     xend = x[4], yend = y[4], colour = age), #neonate
               colour = "gray35", size = 0.8, linejoin = 'mitre', arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))+
  geom_segment(data = trajectory_pcscores_means, aes(x = x[2], y = y[2],  #earlyF
                                                     xend =  x[3], yend = y[3], colour = age), #lateF
               colour = "gray35", size = 0.8, linejoin = 'mitre', arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))
  

#SYMMETRY ANALYSIS ----
#Create links between pairs of symmetric landmarks across line of symmetry to then run symmetry analysis
landpairs <- define.links(mean_shape, ptsize = 4) 

#Analysis of bilateral symmetry to then work on symmetric component alone
minke_sym <- bilat.symmetry(minke_coords,ind = minke_ind, object.sym = TRUE, land.pairs = landpairs) 

#Check ANOVA results - significance value of analysis - if no significant value for "side" means objects are symmetrical
summary(minke_sym)

#Create object to work on symmetric component alone
minke_sym_shape <- minke_sym$symm.shape 

#Mean shape of symmetric component, to use for analysis
mean_shape_sym <- mshape(minke_sym_shape) 

#Plot on mesh to visualize asymmetry
plot(minke_sym, warpgrids = TRUE, mesh = ref_mesh) 
#Save 3D window as html file - 3D widget
PC1min <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

##Perform allometric correction and PCA on symmetry-only shapes
##Allometry regression
minkeAllometry_sym_log <- procD.lm(minke_sym_shape~logCsize, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(minkeAllometry_sym_log) 

#Create residuals array to then save as coordinates for analyses
shape_residuals_sym <- arrayspecs(minkeAllometry_sym_log$residuals,p=dim(minke_sym_shape)[1], k=dim(minke_sym_shape)[2]) 

#New shapes adjusted for allometry with CS to use in analyses
minkeAllometry_sym_residuals <- shape_residuals_sym + array(mean_shape_sym, dim(shape_residuals_sym)) 

#Save mean shape of allometry-adjusted shapes to sue later
mean_shape_sym_res <- mshape(minkeAllometry_sym_residuals)

##Plot shape vs logCS to visualize allometry
#Diagnostic plots to check if model is appropriate - similar to ANOVA tables
init <- par(no.readonly=TRUE) #store initial plot parameters to restore later
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
allometryplot_sym_diagnostics <- plot(minkeAllometry_sym_log,type = "diagnostics",
                                  cex = 1.2, font.main = 2)
par(init)                     #restore initial plot parameters (1 plot showing at a time)


#Regression score of shape vs logCS - regression method with regression score plotting
allometryplot_sym_regscore <- plot(minkeAllometry_sym_log,type = "regression",predictor = logCsize, reg.type = "RegScore",
                               main = "Symmetric shape vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
minkeAllometry_sym_plot <- data.frame(logCS = allometryplot_sym_regscore[["plot.args"]][["x"]], 
                                      RegScores = allometryplot_sym_regscore[["plot.args"]][["y"]])

#Convert data frame to tibble
minkeAllometry_sym_plot <- as_tibble(minkeAllometry_sym_plot)
glimpse(minkeAllometry_sym_plot)
#Add labels and other attributes to tibble as columns
minkeAllometry_sym_plot <- minkeAllometry_sym_plot %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(minkeAllometry_sym_plot)

#Nice plot with specimens colored by age
ggplot(minkeAllometry_sym_plot, aes(x = logCS, y = RegScores, label = individuals, colour = age))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Symmetric shape vs logCS")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

##Add regression line with confidence intervals
#Create object to use for linear model
reg_scores_sym <- allometryplot_sym_regscore[["RegScore"]] 

#Linear model for line
reg_line_sym <- (lm(reg_scores_sym~logCsize))

#Add confidence intervals
#Create data for confidence intervals
x_vals <- seq(min(logCsize), max(logCsize), length = 12)   #use min and max of x values (logCS) as limits and use number of specimens as length of sequence
newX <- expand.grid(logCsize = x_vals)                     #warp x_vals on values of x axis (logCS)
newY <- predict(reg_line_sym, newdata = data.frame(x = newX), interval="confidence",
                level = 0.95)                                      #predict the y values based on the x sequence

#Make data frame of data for confidence intervals
conf_intervals_sym <- data.frame(newX, newY)
#Rename columns to match main plot tibble varibales for x and y
conf_intervals_sym <- rename(conf_intervals_sym, logCS = logCsize, RegScores = fit)
conf_intervals_sym

#Convert data frame to tibble
conf_intervals_sym <- as_tibble(conf_intervals_sym)
glimpse(conf_intervals_sym)
#Add labels and other attributes to tibble as columns to match main plot tibble
conf_intervals_sym <- conf_intervals_sym %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(conf_intervals_sym)

#Nice plot with regression line and confidence intervals - do not color by age or it will mess it up
ggplot(minkeAllometry_sym_plot, aes(x = logCS, y = RegScores, label = individuals))+
  geom_point(size = 3, colour = "chartreuse4")+   #colour all points the same
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Symmetric shape vs logCS")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_smooth(data = conf_intervals_sym, aes(ymin = lwr, ymax = upr), stat = 'identity',     #confidence intervals and reg line
              colour = "darkseagreen3", fill = 'gainsboro')+                             #line colour and interval fill
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1)     #position of tables relative to point (proximity and distance)

##New PCA plot with data corrected for allometry and symmetry
PCA_sym_residuals <- gm.prcomp(minkeAllometry_sym_residuals) 

#List of PC components and proportion of variations
PCA_sym_residuals

#Save PC scores as object to use later
minke_pcscores_sym_res <- PCA_sym_residuals$x

#Save shapes of extremes for axes used in plot
PC1min_sym_res <- PCA_sym_residuals[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_sym_res <- PCA_sym_residuals[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_sym_res <- PCA_sym_residuals[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_sym_res <- PCA_sym_residuals[["shapes"]][["shapes.comp2"]][["max"]] 

#Show deformation grids on axis from mean shape, do this for all 4 extremes - "TPS" method
plotRefToTarget(mean_shape_sym_res, PC1min_sym_res, method = "TPS", mag = 1, label = FALSE)  #save image

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(mean_shape_sym_res, PC1min_sym_res, mesh = ref_mesh, method = "surface", mag = 1, label = FALSE)   #save as HTML

##3D windows save
#Save 3D window as html file - 3D widget
PC1min <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

##Make better PCA plot using ggplot
#Read PC scores as tibble
minke_pcscores_sym_res <- as_tibble(minke_pcscores_sym_res)
glimpse(minke_pcscores_sym_res)
#Add labels and other attributes to tibble as columns
minke_pcscores_sym_res <- minke_pcscores_sym_res %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(minke_pcscores_sym_res)

#Nice plot
ggplot(minke_pcscores_sym_res, aes(x = Comp1, y = Comp2, label = individuals, colour = age))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+            #legend and color adjustments
  theme_bw()+
  xlab("PC 1 (47.38%)")+ #copy this from standard PCA plot or from PCA summary
  ylab("PC 2 (23.76%)")+
  ggtitle("PCA symmetry residuals")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position


#MODULARITY TEST AND INTEGRATION TEST ----
#Set all landmarks in one module
modules <- rep('a',16) 

#Put selected landmarks in second module
modules[4:9]<-'b' 
modules[14:16]<-'b' 
modules

#Perform modularity test - compare selected modules to the null assumption of random assignment of partitions (no modularity at all) 
#Is the data modular according to my defined modules?
#Best done on allometric residuals - NOT ON SYMMETRY data, not a biological significant hypothesis and little difference
modularity_test <- modularity.test(minkeAllometry_residuals, modules, CI = T,iter = 999, print.progress = T) 

#Get P value for CR values (same as RV)
summary(modularity_test) 

#Histogram of CR values
plot(modularity_test)

##Make better histogram plot with ggplot
#Save CR of data as object
CRdata <- modularity_test[["CR"]]

#Save all generated CRs as object
CR <- modularity_test[["CR.boot"]]

#Calculate mean of CRs for arrow
mean_CR <- mean(CR)

#Create tibble with modularity analysis values
modularity_plot <- data.frame(CR)
modularity_plot <- as_tibble(modularity_plot )
glimpse(modularity_plot)

#Make simple histogram plot as object to obtain counts per bin with given bin width (try to see if binwidth value appropiate)
plot1 <- ggplot(modularity_plot, aes(CR))+
  geom_histogram(binwidth = 0.01, fill = "azure2", colour = "azure4")
plot1

#Create data frame with plot variables
plot2 <- ggplot_build(plot1)  

#Save only data counts as tibble
plot2_data <- plot2[["data"]][[1]]
plot2_data <- plot2_data[,1:5]
plot2_data <- as_tibble(plot2_data)

#Filter rows to select row with count for CR data and mean CR 
count_CRdata <- plot2_data %>% filter(xmin <= CRdata, xmax >= CRdata)
count_meanCR <- plot2_data %>% filter(xmin <= mean_CR, xmax >= mean_CR)

#Create tibble with x and y columns to build arrows - x = mean position on bin, y = starting position on bin
arrow_plot <- data.frame(x_data = count_CRdata$x, y_data = count_CRdata$y, x_mean = count_meanCR$x, y_mean = count_meanCR$y)
arrow_plot <- as_tibble(arrow_plot)
glimpse(arrow_plot)

#Nice plot  
ggplot(modularity_plot, aes(CR))+
  geom_histogram(binwidth = 0.01, fill = "azure2", colour = "azure4")+
 #add arrow for CR data
  geom_segment(data = arrow_plot, aes(x = x_data, xend = x_data, y = y_data, yend = y_data + 10, colour = "firebrick4"), size = 1.1,
               arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "first", type = "closed"), linejoin = 'mitre')+
 #add arrow for mean CR 
  geom_segment(data = arrow_plot, aes(x = x_mean, xend = x_mean, y = y_mean, yend = y_mean + 10, colour = "black"), size = 1.1,
               arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "first", type = "closed"), linejoin = 'mitre')+
  scale_colour_manual(name = NULL, labels = c("CR data", "mean CR"), 
                      values = c("firebrick4","black"))+            #legend and color adjustments
  theme_minimal()+
  xlab("CR coefficient")+
  ylab("Frequency")+
  ggtitle("Modularity analysis - p-value=0.0003")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

##Perform integration test between the 2 modules (two-block PLS within configuration) 
#How much are the two or more modules integrated with each other?
#Best done on allometric residuals - NOT ON SYMMETRY data, not a biological significant hypothesis and little difference
integration_test <- integration.test(minkeAllometry_residuals, partition.gp = modules, iter = 999, print.progress = T) 

#Get P value and effect size Z - if significant means the modules are integrated!
summary(integration_test) 

#PLS plot of the two modules - two-block PLS with regression line
modulesplot_pls <- plot(integration_test, 
                        pch = 21, col = "darkgoldenrod2", bg = "darkgoldenrod2", cex = 1.2, ) 

#Save plot arguments as objects to use in plots
block1_modules <- modulesplot_pls$plot.args$x
block2_modules <- modulesplot_pls$plot.args$y

#Add labels
text(x = block1_modules, y = block2_modules, labels = minke_ind,
     pos = 3, offset = 0.5, cex = 0.75)   #improve appearance of labels


##Obtain shapes of modules at min and max of Block 1 and Block 2 axes

#Quick way to visualize shape changes on axes with interactive plot with "points" method
picknplot.shape(modulesplot_pls, label = TRUE, mag = 1)

#Better way using shape.predictor and plotReftoTarget - "TPS" and "surface" methods
#Save for block 1 and block 2 as objects from plot
block1_shapes <- modulesplot_pls$A1
block2_shapes <- modulesplot_pls$A2

#Create mean shape for each block
block1_mean_shape <- mshape(modulesplot_pls$A1)
block2_mean_shape <- mshape(modulesplot_pls$A2)

#Find specimen closer to mean for each block, useful to create warp mesh
findMeanSpec(block1_shapes)
findMeanSpec(block2_shapes)

#Create object containing only that specimen coordinates
block1_mean_spec <- block1_shapes[,,6]
block2_mean_spec <- block2_shapes[,,12]

#Import reference meshes for each block to use in surface method - parts of original reference mesh
block1_3D <- read.ply("Data/simpleskull_block1_rostrum.ply") #make sure NO binary encoding (ASCII)
block2_3D <- read.ply("Data/simpleskull_block2_braincase.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates to make sure it has same scale
range(block1_3D$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(block1_mean_spec)
range(block2_3D$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(block2_mean_spec)

#Create warp meshes, to use as reference for visualization of analyses
block1_ref_mesh <- warpRefMesh(mesh = block1_3D, mesh.coord = block1_mean_spec, ref = block1_mean_shape, color = NULL, centered = FALSE) 
block2_ref_mesh <- warpRefMesh(mesh = block2_3D, mesh.coord = block2_mean_spec, ref = block2_mean_shape, color = NULL, centered = FALSE) 

##Use shape predictor to generate shape coordinate at min and max of each block 
block1_shape_pred <- shape.predictor(block1_shapes, #shapes of block
                             x = block1_modules, #values of axis form plot
                             Intercept = FALSE,  
                             method = "PLS",  #since it's the results of a PLS method
                             min = min(block1_modules), #min value of axis
                             max = max(block1_modules)) #max value of axis

#Show deformation grids on axis from mean shape, do this for all 4 extremes - "TPS" method
plotRefToTarget(block1_mean_shape, block1_shape_pred$min, mag = 1, method = "TPS") #save image
plotRefToTarget(block1_mean_shape, block1_shape_pred$max, mag = 1, method = "TPS") #save image

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(block1_mean_shape, block1_shape_pred$min, mesh = block1_ref_mesh, method = "surface", mag = 1)   #save as HTML
plotRefToTarget(block1_mean_shape, block1_shape_pred$max, mesh = block1_ref_mesh, method = "surface", mag = 1)   #save as HTML

##3D windows save
#Save screenshot of 3D window, useful for lateral and dorsal views - use screen snip if it fails
rgl.snapshot(filename = "Output/X.png") 
#Save 3D window as html file - 3D widget
block1_pred <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)

#Repeat for other block
block2_shape_pred <- shape.predictor(block2_shapes, #shapes of block
                                     x = block2_modules, #values of axis form plot
                                     Intercept = FALSE,  
                                     method = "PLS",  #since it's the results of a PLS method
                                     min = min(block2_modules), #min value of axis
                                     max = max(block2_modules)) #max value of axis

#Show deformation grids on axis from mean shape, do this for all 4 extremes - "TPS" method
plotRefToTarget(block2_mean_shape, block2_shape_pred$min, mag = 1, method = "TPS") #save image
plotRefToTarget(block2_mean_shape, block2_shape_pred$max, mag = 1, method = "TPS") #save image

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(block2_mean_shape, block2_shape_pred$min, mesh = block2_ref_mesh, method = "surface", mag = 1)   #save as HTML
plotRefToTarget(block2_mean_shape, block2_shape_pred$max, mesh = block2_ref_mesh, method = "surface", mag = 1)   #save as HTML

##3D windows save
#Save screenshot of 3D window, useful for lateral and dorsal views - use screen snip if it fails
rgl.snapshot(filename = "Output/X.png") 
#Save 3D window as html file - 3D widget
block2_pred <- scene3d()
widget <- rglwidget()
filename <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(rglwidget(), filename)
browseURL(filename)    #from browser save screenshots as PNG (right click on image-save image) and save HTML (right click on white space-save as->WebPage HTML, only)


##Make better PLS plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
modules_pls_plot <- data.frame(block1_modules, block2_modules)
modules_pls_plot

#Convert data frame to tibble
modules_pls_plot <- as_tibble(modules_pls_plot)
glimpse(modules_pls_plot)
#Add labels and other attributes to tibble as columns
modules_pls_plot <- modules_pls_plot %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(modules_pls_plot)

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(modules_pls_plot, aes(x = block1_modules, y = block2_modules, label = individuals, colour = age))+
  #confidence intervals and reg line, before points  
  geom_smooth(method='lm',     #use standard function, values too small    
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+  
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  xlab("PLS1 Block 1: rostrum")+
  ylab("PLS1 Block 1: braincase")+
  ggtitle ("PLS1 plot: Block 1 (rostrum) vs Block 2 (braincase)")+
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5), axis.title = element_text(size = 11))


#PHYLOGENTIC/AGE ANALYSIS ----
#Import trees in Nexus format - branch lenghts needed!!
trees <- "Data/trees_age_2.nex" 
trees2 <- "Data/trees_age.nex"

##Read the trees for analysis
age_trees <- read.nexus(trees) 
plot(age_trees)
View(age_trees)
spec_trees <- read.nexus(trees2)
plot(spec_trees)
View(spec_trees)

#Create mean shape for each group based on classifiers
age_means <- aggregate(two.d.array(minke_coords) ~ minke_age, FUN = mean) #use original shapes, NOT allometric correction

#Set row names that match tree labels
rownames(age_means) <- c("adult","earlyFetus","lateFetus","neonate") 
age_means

#Transform in matrix and exclude columns that do not contain coordinates
age_means <- as.matrix(age_means[,-1]) 
#Transform in 3D array
age_means <- arrayspecs(age_means,16,3) 

#Mean logCS size based on groups
logCsize_age_means <- aggregate(logCsize ~ minke_age, FUN = mean) 

#Set row names that match tree labels
rownames(logCsize_age_means) <- c("adult","earlyFetus","lateFetus","neonate") 

#Eliminate first column with names of classifiers
logCsize_age_means <- select(logCsize_age_means,-"minke_age") 

#Make numeric
logCsize_age_means <- as.matrix(logCsize_age_means) 

##Phylomorphospace
#Simple phylogeny projected on morphospace, no calculations
PCA_phylomor_spec1 <- gm.prcomp(minkeAllometry_residuals, phy = spec_trees$age2)
PCA_phylomor_spec2 <- gm.prcomp(minkeAllometry_residuals, phy = spec_trees$age4)

#View results
summary(PCA_phylomor_spec1)
summary(PCA_phylomor_spec2)

#Plot PCA
plot(PCA_phylomor_spec1, phylo = TRUE, #add phylogeny
     main = "PCA phylomorphospace",  pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 
plot(PCA_phylomor_spec2, phylo = TRUE, #add phylogeny
     main = "PCA phylomorphospace",  pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 

#3D phylo plot
plot(PCA_phylomor_spec1, time.plot = TRUE, pch = 21, col = "deeppink", bg = "deeppink", cex = 2, 
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


##phyloPCA
#Phylogeny is used to calculate principal components of variation - GLS method
PCA_phylo_spec1 <- gm.prcomp(minkeAllometry_residuals, phy = spec_trees$age2, GLS = TRUE)
PCA_phylo_spec2 <- gm.prcomp(minkeAllometry_residuals, phy = spec_trees$age4, GLS = TRUE)

#View results
summary(PCA_phylo_spec1)
summary(PCA_phylo_spec2)

#Plot PCA
plot(PCA_phylo_spec1, phylo = TRUE, #add phylogeny
     main = "phyloPCA", pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 
plot(PCA_phylo_spec2, phylo = TRUE, #add phylogeny
     main = "phyloPCA", pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 



##PaCA: Phylogenetically-aligned PCA
#Projection of phylogenetic signal in the first two components, helps to highlight importance of phylogeny relative to other signals
PCA_PaCA_spec1 <- gm.prcomp(minkeAllometry_residuals, phy = spec_trees$age2, align.to.phy = TRUE)
PCA_PaCA_spec2 <- gm.prcomp(minkeAllometry_residuals, phy = spec_trees$age4, align.to.phy = TRUE)

#View results
summary(PCA_PaCA_spec1)
summary(PCA_PaCA_spec2)

#Plot PCA
plot(PCA_PaCA_spec1, phylo = TRUE, #add phylogeny
     main = "PaCa", pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 
plot(PCA_PaCA_spec2, phylo = TRUE, #add phylogeny
     main = "PaCa", pch = 21, col = "deeppink", bg = "deeppink", cex = 1, font.main = 2) 


#Add code for GGplot graphics




#ALLOMETRY ANALYSIS BY GROUP ----
minkeAllometry_log_age <- procD.lm(minke_coords ~ logCsize * minke_age, iter=999, print.progress = TRUE) 
View(minkeAllometry_log_age)

#Main results of ANOVA analysis of allometry with logCS
summary(minkeAllometry_log_age) 

##Plot shape vs logCS to visualize allometry
#Diagnostic plots to check if model is appropriate - similar to ANOVA tables
init <- par(no.readonly=TRUE) #store initial plot parameters to restore later
par(mfrow = c(2, 2))          #arrange all the 4 plots next to each other
allometryplot_diagnostics_age <- plot(minkeAllometry_log_age,type = "diagnostics",
                                   cex = 1.2, font.main = 2)
par(init)                     #restore initial plot parameters (1 plot showing at a time)

#Regression score of shape vs logCS - regression method with regression score plotting
allometryplot_regscore_age <- plot(minkeAllometry_log_age, type = "regression", predictor = logCsize, reg.type = "RegScore",
                                main = "Shape vs logCS by age",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics
text(x = logCsize, y = allometryplot_regscore_age$RegScore, labels = minke_ind,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

##Add regression line with confidence intervals to plot
#Create object to use for linear model
reg_scores_age <- allometryplot_regscore_age[["RegScore"]] 

#Linear model for line
reg_line_age <- (lm(reg_scores~logCsize))
summary(reg_line_age)

#Draw line on plot
abline(reg_line_age, col = "darkseagreen3", 
       lty = 2, lwd = 2)     #line type (e.g. dashed) and width

#Add confidence intervals
#Create data for confidence intervals
x_vals <- seq(min(logCsize), max(logCsize), length = 12)   #use min and max of x values (logCS) as limits and use number of specimens as length of sequence
newX <- expand.grid(logCsize = x_vals)                     #warp x_vals on values of x axis (logCS)
newY <- predict(reg_line_age, newdata = data.frame(x = newX), interval="confidence",
                level = 0.95)                                      #predict the y values based on the x sequence

#Draw confidence intervals lines on plot
matlines(newX, newY[,2:3],                                 #first column of newY not useful, it is the fit, 2 and 3 are the min and max values
         col = "darkseagreen3", lty=1)                     #line graphics

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
minkeAllometry_plot_age <- data.frame(logCS = allometryplot_regscore_age[["plot.args"]][["x"]], RegScores = allometryplot_regscore_age[["plot.args"]][["y"]])
minkeAllometry_plot_age

#Convert data frame to tibble
minkeAllometry_plot_age <- as_tibble(minkeAllometry_plot_age)
glimpse(minkeAllometry_plot_age)
#Add labels and other attributes to tibble as columns
minkeAllometry_plot_age <- minkeAllometry_plot_age %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(minkeAllometry_plot_age)

##Add regression line with confidence intervals
#Make data frame of data for confidence intervals
conf_intervals <- data.frame(newX, newY)
conf_intervals 
#Rename columns to match main plot tibble varibales for x and y
conf_intervals <- rename(conf_intervals, logCS = logCsize, RegScores = fit)
conf_intervals 

#Convert data frame to tibble
conf_intervals <- as_tibble(conf_intervals)
glimpse(conf_intervals)
#Add labels and other attributes to tibble as columns to match main plot tibble
conf_intervals <- conf_intervals %>% mutate(individuals = minke_age_tibble$specimenID, age = minke_age_tibble$age)
glimpse(conf_intervals)

#Nice plot with specimens colored by age AND regression line with confidence intervals
ggplot(minkeAllometry_plot_age, aes(x = logCS, y = RegScores, label = individuals, colour = age))+
  geom_smooth(data = conf_intervals, aes(ymin = lwr, ymax = upr), stat = 'identity',          #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Shape vs logCS by age")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1)     #position of tables relative to point (proximity and distance) 



