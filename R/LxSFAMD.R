# FAMD -Factor analysis for mixed data- (like a PCA but for categorical and coninuos mixed data) for Lx epidemiology (Lx dataset)
#aim: to find out which variables explain more the variability of our data
# project start: 02.11.2017
#Modified (including infor ob temperature, salinity and pH): 30.05.2018

install.packages("FactoMineR")
library("FactoMineR")
install.packages("factoextra")
library(ggplot2)
library("factoextra")

setwd("C:/Users/sa04pm/Desktop/Laminariocolax/")
Lx<-read.table("Lx/LxFAMD2.txt", header=TRUE,fill= TRUE) #windows SAMS


Lx<-read.table("/Users/Pedro/Lx/LxSFAMD.txt", header=TRUE,fill= TRUE) #Mac

head(Lx, 3)

Lx$Month = factor(Lx$Month) #change (replace) month from continuous to categorical variable

res.famd <- FAMD(Lx, graph = FALSE) #To compute FAMD

print(res.famd)

#Eigenvalues / Variances 

library("factoextra")
eig.val <- get_eigenvalue(res.famd)
head(eig.val)

#The function fviz_eig() or fviz_screeplot() [factoextra package] 
#can be used to draw the scree plot 
#(the percentages of inertia explained by each FAMD dimensions)

fviz_screeplot(res.famd)

#Graph of variables
#All variables
#The function get_mfa_var() [in factoextra] is used to extract the results for variables. 
#By default, this function returns a list containing the coordinates, the cos2 and the contribution of all variables:

var <- get_famd_var(res.famd)
var

#The different components can be accessed as follow:

# Coordinates of variables
head(var$coord)
# Cos2: quality of representation on the factore map
head(var$cos2)
# Contributions to the  dimensions
head(var$contrib)

#The following figure shows the correlation between variables - both quantitative and qualitative variables - 
#and the principal dimensions, as well as, the contribution of variables to the dimensions 1 and 2. 
#The following functions [in the factoextra package] are used

#The red dashed line on the graph above indicates the expected average value, If the contributions were uniform. 

# Plot of variables
fviz_famd_var(res.famd, repel = TRUE)
# Contribution to the first dimension
fviz_contrib(res.famd, "var", axes = 1)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2)

##Quantitative variables
#To extract the results for quantitative variables, type this:

quanti.var <- get_famd_var(res.famd, "quanti.var")
quanti.var 

#In this section, we’ll describe how to visualize quantitative variables. 
#Additionally, we’ll show how to highlight variables according to either 
#i) their quality of representation on the factor map or ii) their contributions to the dimensions.

#The R code below plots quantitative variables. We use repel = TRUE, to avoid text overlapping.

fviz_famd_var(res.famd, "quanti.var", repel = TRUE,
              col.var = "black")

#Briefly, the graph of variables (correlation circle) shows the relationship between variables, 
#the quality of the representation of variables, as well as, the correlation between variables and the dimensions. 

#The most contributing quantitative variables can be highlighted on the scatter plot using the argument col.var = "contrib". 
##This produces a gradient colors, which can be customized using the argument gradient.cols.

fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

#Similarly, you can highlight quantitative variables using their cos2 values 
#representing the quality of representation on the factor map. 
#If a variable is well represented by two dimensions, the sum of the cos2 is closed to one. For some of the items, 
#more than 2 dimensions might be required to perfectly represent the data.

# Color by cos2 values: quality on the factor map
fviz_famd_var(res.famd, "quanti.var", col.var = "cos2",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
              repel = TRUE)

#Graph of qualitative variables

quali.var <- get_famd_var(res.famd, "quali.var")
quali.var 

#To visualize qualitative variables, type this:

fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

#Graph of individuals
##NOT APPLICABLE HERE##

ind <- get_famd_ind(res.famd)
ind

fviz_famd_ind(res.famd, col.ind = "cos2", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

#In the plot above, the qualitative variable categories are shown in black. 
#If you don’t want to show them on the plot, use the argument invisible = "quali.var".

fviz_famd_ind(res.famd, col.ind = "cos2", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE,
              invisible = "quali.var")



fviz_mfa_ind(res.famd, 
             habillage = "exp", # color by groups 
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE # Avoid text overlapping
) 

fviz_mfa_ind(res.famd, 
             habillage = "month", # color by groups 
             palette = c("#99CC00", "#99FF00", "#33FF00", "#66CCFF", "#6699FF", "#0000FF", "#FF9933", "#CC6600", "#993300", "#FF0033", "#FF0066", "#FF0099"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE # Avoid text overlapping
) 


#If you want to color individuals using multiple categorical variables at the same time, use the function fviz_ellipses() [in factoextra] as follow:

fviz_ellipses(res.famd, c("exp", "month"), repel = TRUE)


### Correlation matrix ###
#subsetting numeric variables

newdata <- Lx[c(-1,-2,-3,-8,-10)]

cormat <- round(cor(newdata),2)
head(cormat)


library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)

#The function geom_tile()[ggplot2 package] is used 
#to visualize the correlation matrix:

#ugly graph
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

#better graph

#Get the lower and upper triangles of the correlation matrix
#Note that, a correlation matrix has redundant information. 
#We'll use the functions below to set half of it to NA.

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri

#Melt the correlation data and drop the rows with NA values:
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

#organizing by correlation values

# Reorder the correlation matrix


reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd) # hclust for hierarchical clustering order is used 
  cormat <-cormat[hc$order, hc$order]
}

#Reordered correlation data visualization 
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

##add corr values
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  scale_y_discrete(breaks=c("Sal","Density","Temp","Gall.prev","N.stipe","D.stipe","T.length","D.holdfast","Rep.prev","ph"), 
                   labels=c("Salinity", "Density", "Temperature","Gall prevalence","Stipe number","Stipe diameter","Total length","Holdfast diameter","Reproductive individuals","pH")) +
  scale_x_discrete(breaks=c("Sal","Density","Temp","Gall.prev","N.stipe","D.stipe","T.length","D.holdfast","Rep.prev","ph"), 
                   labels=c("Salinity", "Density", "Temp.","Gall prev.","Stipe N","Stipe D","T. length","Holdfast D","Rep. indiv.","pH")) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
  