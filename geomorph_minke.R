#always do this first!!
library(geomorph) 
library(geiger)
library(dplyr)
library(ggplot)
library(tidyr)
library(tidyverse)
library(ggrepel)

#DATA IMPORT AND PREP ----

#Import data into R, make sure only first columns are not numbers
minke_raw <- read_csv("Data/16minkeonly.csv") 
glimpse(minke_raw)

#Average landmark takes
minke_avg_takes<-minke_raw %>% group_by(specimenID) %>% summarize(across(starts_with("Raw"), list(mean)))
glimpse(minke_avg_takes)

#Create data frame with only numerical values and set row names as specimen ID
minke_avg_takes <- remove_rownames(minke_avg_takes)
has_rownames(minke_avg_takes) #steps for tibbles
minke_avg_takes2 <- data.frame(column_to_rownames(minke_avg_takes, var = "specimenID")) #adding row names to data frame
minke_avg_takes <- minke_avg_takes2 
remove(minke_avg_takes2) #change name and remove extra object

#Save individual names as object
minke_ind<-row.names(minke_avg_takes) 

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
#Make data frame for analysis in geomorph
minke_age <- data.frame(column_to_rownames(minke_age_tibble, var = "specimenID")) #adding row names to the tibble and making it a data frame
minke_age <- factor(minke_age) #transform data to factors

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
plotRefToTarget(mean_shape, PC1min_all, method = "TPS", mag = 1, label = FALSE)

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - "points" method
plotRefToTarget(mean_shape, PC1min_all, method = "points", mag = 1, label = FALSE)

#Show 3D deformation from mean with vectors, do this for all 4 extremes - "vector" method
plotRefToTarget(mean_shape, PC1min_all, method = "vector", mag = 1, label = FALSE)

#Show 3D deformation from mean by warping 3D mesh, do this for all 4 extremes - "surface" method
plotRefToTarget(mean_shape, PC1min_all, mesh = ref_mesh, method = "surface", mag = 1, label = FALSE)

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

#Nice plot with specimens colored by age
ggplot(minkeAllometry_plot, aes(x = logCS, y = RegScores, label = individuals, colour = age))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Adult", "Early Fetus", "Late Fetus", "Neonate"), 
                      values = c("blue4","cyan2","deepskyblue1","dodgerblue3"))+           
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Shape vs logCS")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

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

#Nice plot with regression line and confidence intervals - do not color by age or it will mess it up
ggplot(minkeAllometry_plot, aes(x = logCS, y = RegScores, label = individuals))+
  geom_point(size = 3, colour = "chartreuse4")+   #colour all points the same
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Shape vs logCS")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_smooth(data = conf_intervals, aes(ymin = lwr, ymax = upr), stat = 'identity',     #confidence intervals and reg line
              colour = "darkseagreen3", fill = 'gainsboro')+                             #line colour and interval fill
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

#TWO-BLOCK PLS ----


allom.pls<-two.b.pls(logCSsize,minke.shape,iter = 999) #Two-block PLS of allometry, another way to visualize connection between logCS and shape, main one for analyses
summary(allom.pls) #get P-value of regression
plot(allom.pls,pch=19,col = as.numeric(classifiers)) #Make sure it is ok and save plot with pred line


allom.pls.CS<-two.b.pls(CSsize,minke.shape,iter = 999) #Two-block PLS of allometry, another way to visualize connection between CS and shape
summary(allom.pls.CS) #get P-value of regression
plot(allom.pls.CS,pch=19,col = as.numeric(classifiers)) #Make sure it is ok and save plot with pred line

#PCA ALLOMETRY RESIDUALS ----

PCAplot.all<-gm.prcomp(minke.shape.all, label = minke.ind,groups = classifiers,legend = TRUE) #New PCA plot with data corrected for allometry


PCAplot.all #List of PC components and proportion of variations


summary(PCAplot.all) #Vie) results


plot(PCAplot.all) #View plot
PC1min.all<-PCAplot.all$pc.shapes$PC1min #Save shapes of extremes for axes used in plot
PC1max.all<-PCAplot.all$pc.shapes$PC1max #Save shapes of extremes for axes used in plot
PC2min.all<-PCAplot.all$pc.shapes$PC2min #Save shapes of extremes for axes used in plot
PC2max.all<-PCAplot.all$pc.shapes$PC2max #Save shapes of extremes for axes used in plot
PC1min.all.plot<-plotRefToTarget(minke.shape.mean,PC1min.all,method = "TPS",mag = 1, links = links,label = FALSE) #Obtain grids for deformation on axis from mean shape, do this for all 4 extremes
PC1min.all.plot.3D<-plotRefToTarget(minke.shape.mean,PC1min.all,method = "points",mag = 1, links = links,label = FALSE) #Show 3D deformation from mean, do this for all 4 extremes
rgl.snapshot(filename = "X.png") #save screenshot of 3D deformation plot, useful for lateral and dorsal views
#Better use screen snip this fails!
PC1min.all.plot.mesh<-plotRefToTarget(minke.shape.mean,PC1min.all,method = "surface",mesh=minke.skull3D, mag = 1, label = FALSE) #Show mesh3D that represents extreme of axis, do this for all 4 extremes
writePLY("PC1min.all.mesh.ply") #save as PLY file





#ANOVA OF GROUP DIFFERENCES  ----
dataframe.symm<-geomorph.data.frame(minke.shape,gp=classifiers) #Create dataframe to operate more easily
age.anova<-procD.lm(minke.shape~gp,iter=999,data=dataframe.symm) #Conduct ANOVA to test differences between groups
summary(age.anova) #Results and significance
PCAplot.age.anova<-plotTangentSpace(minke.shape,groups = dataframe.symm$gp,label=minke.ind) #PCA of residuals of ANOVA, highlights groups





#TRAJECTORY ANALYSIS ----
trajectory.age<- trajectory.analysis(minkeAllometry.log,groups = classifiers,traj.pts=2,pca = TRUE,print.progress = TRUE) #Main analysis, use obj from procD.lm
summary(trajectory.age,show.trajectories = TRUE) #View results
trajectory.plot<-plot(trajectory.age) #Plot results
add.trajectories(trajectory.plot,traj.pch = 21,traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, start.bg = 3, end.bg = 2) #Add line between groups

#SYMMETRY ANALYSIS ----

#Create links between landmarks across line of symmetry to then run symmetry analysis
minke.landpairs<-define.links(mean.shape, ptsize = 2) 
#Analysis of bilateral symmetry to then work on symmetric component alone
minke.sym<-bilat.symmetry(minke.allcoords,ind=minke.ind, object.sym = TRUE,land.pairs = minke.landpairs) 
#Check ANOVA results
minke.sym 
#Create object to work on symmetric component alone
minke.shape<-minke.sym$symm.shape 
#Mean shape of symmetric component, to use for analysis
minke.shape.mean<-mshape(minke.shape) 
#plot specimens to make sure it works
plotAllSpecimens(minke.shape,links = links) 

#MODULARITY TEST ----
land.gps<-rep('a',16) #Set all landmarks in one module
land.gps[4:9]<-'b' #Put selected landmarks in second module
land.gps[14:16]<-'b' #Put selected landmarks in second module
MT <- modularity.test(minke.shape.all,land.gps,CI=F,iter=999, print.progress = T) #Perform modularity test
summary(MT) #Get P value for CR values (same as RV)
plot(MT) #Histogram of CR values
IT <- integration.test(minke.shape.all, partition.gp=land.gps, iter=999, print.progress =T) #Perform integration test between the 2 modules (PLS)
summary(IT) #Get P value and PLS results
plot(IT) #PLS plot of the modules

#GROUP MEANS FOR PHYLOGENTIC/AGE ANALYSIS ----
trees<-"trees_age.nex" #Add trees in Nexus format
agetrees<-read.nexus(trees) #Read the trees for analysis
means <- aggregate(two.d.array(minke.shape.all)~ classifiers, FUN=mean) #Create mean shape for each group based on classifiers
rownames(means) <- c("adult","earlyFetus","lateFetus","neonate") #Set row names that match tree labels
means.age <- as.matrix(means[,-(1)]) #Transform in matrix and exclude columns that do not contain coordinates
means.age<-arrayspecs(means.age,16,3) #Transform in 3D array as whole data
means.CS<-aggregate(logCSsize~ classifiers, FUN=mean) #Mean CS size based on groups
rownames(means.CS) <- c("adult","earlyFetus","lateFetus","neonate") #Set row names that match tree labels
means.CS.age<-select(means.CS,-"classifiers") #Eliminate first column with names of classifiers
means.CS.age<-as.matrix(means.CS.age) #Make numeric