##GGPLOT CODE##

library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(borealis)
library(ggthemes)
library(gridExtra)

#Make palette with ggthemes - color or shapes
View(ggthemes_data)

#Chose palette from viewer and send to console, then copy in script
mypalette <- ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Blue"]][["value"]]

#Image entire palette, number of values given by viewer
image(1:20, 1, as.matrix(1:20),col = mypalette, xlab = "Blue (sequential)", ylab = "", yaxt = "n")

#Choose colors based on number of factors/groupings and image to check
mypalette <- as.matrix(mypalette)
project_palette <- c(mypalette[3,], mypalette[9,], mypalette[15,], mypalette[20,])
image(1:4, 1, as.matrix(1:4), col = project_palette, xlab = "Project palette", ylab = "", yaxt = "n")

#Create shape palette
shapes <- c(21, 22, 24) #circle, square and triangle, fillable and with border

#Convert data frame to tibble
pcscores_all_size <- as_tibble(pcscores_all_size)
glimpse(pcscores_all_size)
#Add labels and other attributes to tibble as columns
pcscores_all_size <- pcscores_all_size %>% mutate(individuals = classifiers$specimenID, age = classifiers$age)
#Check tibble is ordered as initial data
glimpse(pcscores_all_size)

#PCA/regression plot with geom_points and scale color ----
ggplot(pcscores_all, aes(x = Comp1, y = Comp2, label = individuals, colour = age))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+            #legend and color adjustments
  theme_bw()+   #theme_light nd theme_minimal also ok
  xlab("PC 1 (48.89%)")+ #copy this from standard PCA plot
  ylab("PC 2 (22.17%)")+
  ggtitle("PCA all data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position

#PCA with borealis gg.shape.space() and ggplot ----
#borealis plot
PCA_all_ggplot <- gg.shape.space(PCA_all, #gm.prcomp object
                                 group = pcscores_all$age, convex.hulls = TRUE,  #divide by group and draw convex hulls
                                 pt.size = 3, label.groups = F,  #point size and delete group labels
                                 color = project_palette, include.legend = F,  #colour as project palette but do not add legend
                                 backtransform.examples = F) #avoid shapes on plot, too crowded

#Nice plot with extra elements from ggplot
PCA_all_ggplot +
  #add correct legend for groups - ignore error message after code runs
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+ 
  #add labels for points  
  geom_text_repel(aes(x = Comp1, y = Comp2, label = individuals),      #aes() from tibble
                  data = pcscores_all, colour = "black", size = 3.5)+  #text size and color for labels
  #change theme
  theme_bw()+
  #add title and change its appearance 
  ggtitle("PCA all data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 


#Regression plot with geom_points and geom_abline with intercept and slope from external model ----

#Calculate regression for each component
reg_PC1all_size <- lm(Comp1 ~ size, data = pcscores_all_size)

#Create data frame with line parameters from regression
#Allows to show a line on PC plot with specimens colored IF col and other graphics OUTSIDE of aes()!!!
regline_PCAall_size <- data.frame(int1 = reg_PC1all_size[["coefficients"]][["(Intercept)"]], 
                                  slope1 = reg_PC1all_size[["coefficients"]][["size"]], #PC1 values
                                  int2 = reg_PC2all_size[["coefficients"]][["(Intercept)"]], 
                                  slope2 = reg_PC2all_size[["coefficients"]][["size"]]) #PC2 values

#Plot
ggplot(pcscores_all_size, aes(x = size, y = Comp1, label = individuals, colour = age))+
  #line on plot
  geom_abline(data = regline_PCAall_size, aes(intercept = int1, slope = slope1), colour = "darkblue", size = 0.8, linetype = "dashed", show.legend = FALSE)+
  geom_point(size = 3)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+           
  theme_classic(base_size = 12)+
  xlab("logCS")+
  ylab("PC 1 (48.89%)")+
  ggtitle("PC1 vs logCS - p-value = 0.007**")+  #copy from summary linear model
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))+
  geom_text_repel(colour = "black", size = 3.5)

#Regression plot with geom_points and geom_smooth with line and confidence intervals from external model ----

##Add regression line with confidence intervals to plot
#Create object to use for linear model
allometry_regscores <- allometry_plot_regscore[["RegScore"]] 

#Linear model for line
allometry_regline <- (lm(allometry_regscores ~ logCsize))

#Add confidence intervals
#Create data for confidence intervals
allometry_xvals <- seq(min(logCsize), max(logCsize), length = 12)   #use min and max of x values (logCS) as limits and use number of specimens as length of sequence
allometry_newX <- expand.grid(logCsize = allometry_xvals)                     #warp x_vals on values of x axis (logCS)
allometry_newY <- predict(allometry_regline, newdata = data.frame(x = allometry_newX), interval="confidence",
                          level = 0.95) 
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

#Plot
ggplot(allometry_plot_tibble, aes(x = logCS, y = RegScores, label = individuals, colour = age))+
  geom_smooth(data = allometry_conf_intervals, aes(ymin = lwr, ymax = upr), stat = 'identity',          #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+           
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Shape vs logCS - p-value = 0.004**")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1)     #position of tables relative to point (proximity and distance) 


#Regression plot with geom_point and geom_smooth with line and confidence intervals using standard function  ----

#Plot
ggplot(allometry_plot_CAC_tibble, aes(x = CAC, y = RSC1, colour = age))+
  #difficult to do with external model for residuals (too small)
  geom_smooth(method = 'lm',        #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.5)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+           
  theme_classic(base_size = 12)+
  ggtitle ("Residual shape component (RSC) 1 vs CAC")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Trajectory analysis plot with geom_points and geom_segment for arrows with coordinates from external data frame ----
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

#Order tibble by avariable e.g. age
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

#Order tibble by avariable e.g. age
#Make factor for variable
group_trajectory_pcscores_means$age <- factor(group_trajectory_pcscores_means$age, 
                                              levels = c("earlyFetus", "lateFetus", "neonate", "adult")) #use the original factor to copy the list of levels
#Order
group_trajectory_pcscores_means <- group_trajectory_pcscores_means[order(group_trajectory_pcscores_means$age),]
glimpse(group_trajectory_pcscores_means)

#Rename columns so they are easier to use for plot
group_trajectory_pcscores_means <- group_trajectory_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
group_trajectory_pcscores_means

#Plot
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
                      values = project_palette)+            #legend and color adjustments
  theme_bw()+
  xlab("PC 1 (60.17%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (30.92%)")+
  ggtitle("Trajectories by age")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Histogram plot using geom_histogram with 2 arrows for data and mean values ----

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
  ggtitle("Modularity analysis - p-value <0.05***")+  #copy data from standard plot
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))




#Combine mutiple plots in 1 - 2 methods ----

#Multiple plots together specifying how many per row or col (ggfortify) - not great, need shared legend
new('ggmultiplot', plots = list(p1, p2), ncol = 1, nrow = 1)  #save the ggplot plots as objects to be listed

#Plot together with shared legend
#Create plot to extract legend
plot_legend <- ggplot(proj_Neoph_tibble, aes(x = axis1, y = axis2, colour = specimens))+
  geom_point(size = 3)+
  scale_colour_manual(name = "Specimens", labels = c("1-Delphinapterus_554", "2-Delphinapterus_555", "3-Delphinapterus_556",
                                                     "4-Megaptera_MVZ", "5-Megaptera_SDNHM", "6-Neophocoena_1889", "7-Physeter_4.8"), #to be ordered as they appear in tibble
                      values = scale)+ 
  guides(colour = guide_legend(override.aes = list(shape = NA)))+
  theme_bw()+
  theme(legend.title.align = 0.5, legend.title = element_text(size = 9), 
        legend.direction = "vertical", legend.position = "right", legend.text = element_text(size = 8))

plot_legend

#Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

shared_legend <- extract_legend(plot_legend)

grid.arrange(arrangeGrob(plot_Meg, plot_Neoph, plot_noNeoph, ncol = 3),
             shared_legend, nrow = 2, heights = c(10, 7))
#Plot with hulls ----
#Create hulls
hulls <- proj_all_tibble_hulls %>%
  group_by(number) %>%
  slice(chull(axis1, axis2)) %>%
  rename(x = axis1, y = axis2)
View(hulls)

#Plot with hulls
ggplot(proj_all_tibble, aes(x = axis1, y = axis2, label = number, colour = specimens,  fill = specimens))+
  geom_point(size = 3)+
  geom_polygon(data = hulls, aes(x = x, y = y), alpha = .5, show.legend = FALSE)+
  theme_bw()+   #theme_light nd theme_minimal also ok
  xlab("I (35.2%-26.3%)")+ #copy this from values tibble
  ylab("II (27.9%-23.3%)")+
  ggtitle("GPSA all tests")+
  guides(alpha = guide_legend(override.aes = list(shape = NA)))+
  geom_text_repel(colour = "black", size = 5, max.overlaps = 20, show.legend = FALSE)+
  theme(legend.title.align = 0.5, legend.direction = "horizontal", legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13))  #title font and position
###Extra code ggplot ----
#Code to remove dots from legend
guides(fill = guide_legend(override.aes = list(shape = NA)))+  #remove annoying dots in the colour legend
  #Code to decide which polygons are sued for each group if shape = present
  scale_shape_manual(values=c(17, 18, 15)) 