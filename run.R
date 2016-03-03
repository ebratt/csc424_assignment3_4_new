##########################
# PROBLEM 1              #
##########################
# setup
# clear the environment
rm(list=ls())

DATA_DIR <- './data'
IMAGES_DIR <- './images'
OUTPUT_DIR <- './output'

make_dir <- function(d) {
    if (file.exists(d)) unlink(d, recursive=TRUE, force=TRUE)
    dir.create(d)
}
lapply(c(IMAGES_DIR, OUTPUT_DIR),make_dir)


## function that concatenates strings (useful for directory paths)
concat <- function(x1,x2) {
    result <- paste(x1,x2,sep="")
    return(result)
}

## function that checks to see if a package is installed and,if not,installs it
## portions of this code came from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
load_package <- function(x) {
    if (x %in% rownames(installed.packages())) { 
        print(concat("package already installed: ", x))
    }
    else { 
        install.packages(x) 
    }
    library(x, character.only=TRUE)
}

# get the data
data <- read.table(concat(DATA_DIR, "/beetle.txt"))  # read text file 
# how many are n/a?
sum(is.na(data))
head(which(is.na(data)))
# how many are NULL?
sum(is.null(data))
# how many are blank?
length(which(data == ""))
str(data)
summary(data)
load_package("psych")
describe(data)

# data dictionary
# V1: flea beetle species; (a) Halticus oleracea (b) Halticus cardourum
# V2: thorax length (THORAX) in microns
# V3: elytra length (ELYTRA) in microns
# V4: length of second antennal joint (AJ2) in microns
# V5: length of third antennal joint (AJ3) in microns

# relabel the variables
colnames(data) <- c("species", "thorax", "elytra", "aj2", "aj3")
# label the classifier
data$species_name <- ifelse(data$species == "a", "oleracea","cardourum")
head(data)

# Look at the profile plot
makeProfilePlot <- function(mylist,names)
{
    # from http://little-book-of-r-for-multivariate-analysis.readthedocs.org/en/latest/src/multivariateanalysis.html
    require(RColorBrewer)
    # find out how many variables we want to include
    numvariables <- length(mylist)
    # choose 'numvariables' random colours
    colours <- brewer.pal(numvariables,"Set1")
    # find out the minimum and maximum values of the variables:
    mymin <- 1e+20
    mymax <- 1e-20
    for (i in 1:numvariables)
    {
        vectori <- mylist[[i]]
        mini <- min(vectori)
        maxi <- max(vectori)
        if (mini < mymin) { mymin <- mini }
        if (maxi > mymax) { mymax <- maxi }
    }
    # plot the variables
    for (i in 1:numvariables)
    {
        vectori <- mylist[[i]]
        namei <- names[i]
        colouri <- colours[i]
        if (i == 1) { plot(vectori,col=colouri,type="l",ylim=c(mymin,mymax)) }
        else         { points(vectori, col=colouri,type="l")                                     }
        lastxval <- length(vectori)
        lastyval <- vectori[length(vectori)]
        text((lastxval-10),(lastyval),namei,col="black",cex=0.6)
    }
}
load_package("RColorBrewer")
names <- colnames(data[,2:5])
mylist <- as.list(data[,2:5])
makeProfilePlot(mylist,names)

# look at means and sd's
sapply(data[,2:5],mean)
sapply(data[,2:5],sd)
# they appear to be consistent when not grouped

# what about grouped by species?
printMeanAndSdByGroup <- function(variables,groupvariable)
{
    # from http://little-book-of-r-for-multivariate-analysis.readthedocs.org/en/latest/src/multivariateanalysis.html
    # find the names of the variables
    variablenames <- c(names(groupvariable),names(as.data.frame(variables)))
    # within each group, find the mean of each variable
    groupvariable <- groupvariable[,1] # ensures groupvariable is not a list
    means <- aggregate(as.matrix(variables) ~ groupvariable, FUN = mean)
    names(means) <- variablenames
    print(paste("Means:"))
    print(means)
    # within each group, find the standard deviation of each variable:
    sds <- aggregate(as.matrix(variables) ~ groupvariable, FUN = sd)
    names(sds) <- variablenames
    print(paste("Standard deviations:"))
    print(sds)
    # within each group, find the number of samples:
    samplesizes <- aggregate(as.matrix(variables) ~ groupvariable, FUN = length)
    names(samplesizes) <- variablenames
    print(paste("Sample sizes:"))
    print(samplesizes)
}
printMeanAndSdByGroup(data[2:5],data[1])

# look at within group variances
calcWithinGroupsVariance <- function(variable,groupvariable)
{
    # from http://little-book-of-r-for-multivariate-analysis.readthedocs.org/en/latest/src/multivariateanalysis.html
    # find out how many values the group variable can take
    groupvariable2 <- as.factor(groupvariable[[1]])
    levels <- levels(groupvariable2)
    numlevels <- length(levels)
    # get the mean and standard deviation for each group:
    numtotal <- 0
    denomtotal <- 0
    for (i in 1:numlevels)
    {
        leveli <- levels[i]
        levelidata <- variable[groupvariable==leveli,]
        levelilength <- length(levelidata)
        # get the standard deviation for group i:
        sdi <- sd(levelidata)
        numi <- (levelilength - 1)*(sdi * sdi)
        denomi <- levelilength
        numtotal <- numtotal + numi
        denomtotal <- denomtotal + denomi
    }
    # calculate the within-groups variance
    Vw <- numtotal / (denomtotal - numlevels)
    return(Vw)
}

# calculate between group variances
calcBetweenGroupsVariance <- function(variable,groupvariable)
{
    # from http://little-book-of-r-for-multivariate-analysis.readthedocs.org/en/latest/src/multivariateanalysis.html
    # find out how many values the group variable can take
    groupvariable2 <- as.factor(groupvariable[[1]])
    levels <- levels(groupvariable2)
    numlevels <- length(levels)
    # calculate the overall grand mean:
    grandmean <- mean(variable[[1]])
    # get the mean and standard deviation for each group:
    numtotal <- 0
    denomtotal <- 0
    for (i in 1:numlevels)
    {
        leveli <- levels[i]
        levelidata <- variable[groupvariable==leveli,]
        levelilength <- length(levelidata)
        # get the mean and standard deviation for group i:
        meani <- mean(levelidata)
        sdi <- sd(levelidata)
        numi <- levelilength * ((meani - grandmean)^2)
        denomi <- levelilength
        numtotal <- numtotal + numi
        denomtotal <- denomtotal + denomi
    }
    # calculate the between-groups variance
    Vb <- numtotal / (numlevels - 1)
    Vb <- Vb[[1]]
    return(Vb)
}

# separation is the between-group variance divided by the within-group variance
calcSeparations <- function(variables,groupvariable)
{
    # find out how many variables we have
    variables <- as.data.frame(variables)
    numvariables <- length(variables)
    # find the variable names
    variablenames <- colnames(variables)
    # calculate the separation for each variable
    for (i in 1:numvariables)
    {
        variablei <- variables[i]
        variablename <- variablenames[i]
        Vw <- calcWithinGroupsVariance(variablei, groupvariable)
        Vb <- calcBetweenGroupsVariance(variablei, groupvariable)
        sep <- Vb/Vw
        print(paste("variable",variablename,"Vw=",Vw,"Vb=",Vb,"separation=",sep))
    }
}
calcSeparations(data[2:5],data[1])

# find the most highly-correlated variables
mosthighlycorrelated <- function(mydataframe,numtoreport)
{
    # find the correlations
    cormatrix <- cor(mydataframe)
    # set the correlations on the diagonal or lower triangle to zero,
    # so they will not be reported as the highest ones:
    diag(cormatrix) <- 0
    cormatrix[lower.tri(cormatrix)] <- 0
    # flatten the matrix into a dataframe for easy sorting
    fm <- as.data.frame(as.table(cormatrix))
    # assign human-friendly names
    names(fm) <- c("First.Variable", "Second.Variable","Correlation")
    # sort and print the top n correlations
    head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
}
mosthighlycorrelated(data[2:5], 4)
cor(data[,2:5])

# scatterplot matrix - no grouping
load_package("car")
png(concat(IMAGES_DIR,'/beetle scatterplot 1.png'), 
    width = 1024, height = 1024)
scatterplotMatrix(data[,2:5], 
                  diagonal="density",
                  main="Scatterplot Matrix for Flea Beetle Species")
dev.off()

# grouped by species
png(concat(IMAGES_DIR,'/beetle scatterplot 2.png'), 
    width = 1024, height = 1024)
scatterplotMatrix(~thorax + elytra + aj2 + aj3|species_name, 
                  data=data, 
                  diagonal="density",
                  main="Scatterplot Matrix for Flea Beetle Species - Grouped by Species")
dev.off()

# color-coded
png(concat(IMAGES_DIR,'/beetle scatterplot 3.png'), 
    width = 1024, height = 1024)
pairs(data[c("thorax", "elytra", "aj2", "aj3")], 
      main="Scatterplot Matrix for Flea Beetle Species - Colored by Species", 
      pch=22, 
      bg=c("red", "blue")[unclass(data$species)])
dev.off()

# create the covariance matrix
cm <- cov(data[2:5])
cov_eigs <- eigen(cm)
cov_eigs$values
cov_eigs$vectors
png(concat(IMAGES_DIR,'/beetle covarience eigenvalues.png'), 
    width = 1024, height = 1024)
plot(cov_eigs$values, type="b")
dev.off()
# create the correlation matrix
cr <- cor(data[2:5])
cor_eigs <- eigen(cr)
cor_eigs$values
cor_eigs$vectors
png(concat(IMAGES_DIR,'/beetle correlation eigenvalues.png'), 
    width = 1024, height = 1024)
plot(cor_eigs$values, type="b")
dev.off()

# Linear Discriminant Analysis
load_package("MASS")
(lda <- lda(species_name ~ ., data=data[,2:6]))

############################
# answer to i.             #
############################
# Box's M: It performs the Box's M-test for homogeneity of covariance matrices 
# obtained from multivariate normal data according to one classification factor. 
# The test is based on the chi-square approximation.
load_package("biotools")
boxM(data[,2:5], grouping=data$species_name)
# Box's M-test for Homogeneity of Covariance Matrices
# 
# data:  data[, -1]
# Chi-Sq (approx.) = 8.4516, df = 10, p-value = 0.5848

############################
# answer to ii.            #
############################
ld_coefs <- lda$scaling
write.table(round(ld_coefs,3), 
            file=concat(OUTPUT_DIR,'/beetle lda coefficients.csv'), sep=",")
# Z = -0.093 * thorax + 0.039 * elytra + 0.024 * aj2 + 0.037 * aj3
# scaled so that their mean value is zero and variance is one

############################
# answer to iii.           #
############################
# make predictions
plda <- predict(lda)
png(concat(IMAGES_DIR,'/beetle lda histograms.png'), 
    width = 1024, height = 1024)
ldahist(data=plda$x[,1], g=data$species_name)
dev.off()
new.beetle <- data.frame(thorax=184,
                         elytra=275,
                         aj2=143,
                         aj3=192)
predict(lda, newdata = new.beetle)
# we would classify this new specimen as oleracea

############################
# answer to iv.            #
############################
# Confusion Matrix:
confusion_matrix <- table(plda$class, data$species_name)
write.table(confusion_matrix, 
            file=concat(OUTPUT_DIR,'/beetle confusion matrix.csv'), sep=",")
# estimate the percentage of beetles that will be mis-classified
round(1 - diag(prop.table(confusion_matrix)), 4)
# total percent incorrect
round(1 - sum(diag(prop.table(confusion_matrix))), 4)
# answer to v.: 
# cross-validation with leave-one-out
lda2 <- lda(species_name ~ ., data=data[,2:6], CV=TRUE)
# Confusion Matrix:
confusion_matrix2 <- table(lda2$class, data$species_name)
write.table(confusion_matrix2, 
            file=concat(OUTPUT_DIR,'/beetle confusion matrix leave one out.csv'), sep=",")
# estimate the percentage of beetles that will be mis-classified
round(1 - diag(prop.table(confusion_matrix2)), 4)
# total percent incorrect
round(1 - sum(diag(prop.table(confusion_matrix2))), 4)

# Just for fun, plot the Partitions
load_package("klaR")
png(concat(IMAGES_DIR,'/beetle projected lda partitionings.png'), 
    width = 1024, height = 1024)
partimat(species ~ thorax + elytra + aj2 + aj3,
         data=data[,1:5],
         method="lda",
         main="Flea Beetle Data Partitioned by Species")
dev.off()


############################
# PROBLEM 3                #
############################
# setup
# clear the environment
rm(list=ls())

DATA_DIR <- './data'
IMAGES_DIR <- './images'
OUTPUT_DIR <- './output'

make_dir <- function(d) {
    if (file.exists(d)) unlink(d, recursive=TRUE, force=TRUE)
    dir.create(d)
}
lapply(c(IMAGES_DIR, OUTPUT_DIR),make_dir)


## function that concatenates strings (useful for directory paths)
concat <- function(x1,x2) {
    result <- paste(x1,x2,sep="")
    return(result)
}

## function that checks to see if a package is installed and,if not,installs it
## portions of this code came from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
load_package <- function(x) {
    if (x %in% rownames(installed.packages())) { 
        print(concat("package already installed: ", x))
    }
    else { 
        install.packages(x) 
    }
    library(x, character.only=TRUE)
}

# get the data
load_package("foreign")
data <- read.spss(concat(DATA_DIR, "/faculty.sav"), 
                  use.value.labels=TRUE,
                  to.data.frame = TRUE,
                  use.missings=FALSE)
keeps <- c("facrank", 
           "item13",
           "item14",
           "item15",
           "item16",
           "item17",
           "item18",
           "item19",
           "item20",
           "item21",
           "item22",
           "item23",
           "item24")
data <- data[keeps]

# Linear Discriminant Analysis with Cross-Validation
load_package("MASS")
x <- data[, 2:13]
data[, 2:13] <- scale(x, 
                      center=FALSE, 
                      scale=apply(x,2, sd, na.rm=TRUE))
head(data)

# Performa Linear Descriminant Analysis
lda <- lda(formula=data$facrank ~ ., data=data)
ld_coefs <- lda$scaling
write.table(round(ld_coefs,3), file=concat(OUTPUT_DIR,'/faculty lda coefficients.csv'), sep=",")
# Box's M: It performs the Box's M-test for homogeneity of covariance matrices 
# obtained from multivariate normal data according to one classification factor. 
# The test is based on the chi-square approximation.
load_package("biotools")
boxM(data[,2:13], grouping=data[,1])
# Box's M test suggests that the variance-covariance matrices are not
# homegeneous:
# Box's M-test for Homogeneity of Covariance Matrices
# 
# data:  data[, 2:13]
# Chi-Sq (approx.) = 572.96, df = 234, p-value < 2.2e-16

# make predictions
plda <- predict(lda)

# do the histograms overlap?
png(concat(IMAGES_DIR,'/faculty ld1 histogram.png'), 
    width = 1024, height = 1024)
ldahist(data=plda$x[,1], g=data$facrank)
dev.off()
png(concat(IMAGES_DIR,'/faculty ld2 histogram.png'), 
    width = 1024, height = 1024)
ldahist(data=plda$x[,2], g=data$facrank)
dev.off()

# Confusion Matrix:
confusion_matrix <- table(plda$class, data$facrank)
write.table(confusion_matrix, 
            file=concat(OUTPUT_DIR,'/faculty confusion matrix.csv'), sep=",")
# estimate the percentage of faculty rankings that will be mis-classified
round(1 - diag(prop.table(confusion_matrix)), 4)
# total percent incorrect
round(1 - sum(diag(prop.table(confusion_matrix))), 4)
# cross-validation with leave-one-out
lda2 <- lda(data$facrank ~ ., data=data[,2:13], CV=TRUE)
# Confusion Matrix:
confusion_matrix2 <- table(lda2$class, data$facrank)
write.table(confusion_matrix2, 
            file=concat(OUTPUT_DIR,'/faculty confusion matrix leave one out.csv'), sep=",")
# estimate the percentage of beetles that will be mis-classified
round(1 - diag(prop.table(confusion_matrix2)), 4)
# total percent incorrect
round(1 - sum(diag(prop.table(confusion_matrix2))), 4)

# plot the LDA projection
prop.lda = lda$svd^2/sum(lda$svd^2)
lda_data <- data.frame(facrank = data[,"facrank"], lda = plda$x)
lda_data <- lda_data[,1:3] # drop unnecessary LD's
load_package("ggplot2")
load_package("scales")
png(concat(IMAGES_DIR,'/faculty lda plot.png'), 
    width = 1024, height = 1024)
ggplot(lda_data) + 
    geom_point(aes(lda.LD1, lda.LD2, col = facrank, shape = facrank), size = 2.5) + 
    labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
         y = paste("LD2 (", percent(prop.lda[2]), ")", sep="")) +
    ggtitle("LDA Projection of Faculty Data") + 
    theme(plot.title = element_text(lineheight=.8, face="bold"))
dev.off()

############################
# answer to a.             #
############################
scaled_data <- data[,2:13]
# Determine number of clusters based on weighted sum of squares
wss <- (nrow(scaled_data)-1)*sum(apply(scaled_data,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(scaled_data, 
                                     centers=i)$withinss)
png(concat(IMAGES_DIR,'/faculty wss to pick number of clusters.png'), 
    width = 1024, height = 1024)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()
# K-Means Cluster Analysis for k=5 on original data
fit5 <- kmeans(scaled_data, 5) # 5 cluster solution

# cluster plot with ellipses
load_package("cluster")
png(concat(IMAGES_DIR,'/faculty cluster plot1.png'), 
    width = 1024, height = 1024)
clusplot(scaled_data, 
         fit5$cluster, 
         color=TRUE, 
         shade=TRUE, 
         labels=4, 
         lines=0,
         main="Cluster Plot of Faculty Data with k=5")
dev.off()

# Centroid Plot against 1st 2 discriminant functions
load_package("fpc")
png(concat(IMAGES_DIR,'/faculty cluster plot2.png'), 
    width = 1024, height = 1024)
plotcluster(scaled_data, 
            method="dc",
            fit5$cluster,
            main="Cluster Plot of Faculty Data with k=5")
dev.off()



############################
# answer to b.             #
############################
d <- dist(scaled_data, method="euclidean")
# Option "ward.D2" implements Ward's (1963) clustering criterion 
# (Murtagh and Legendre 2014). With the latter, the dissimilarities 
# are squared before cluster updating.
fit <- hclust(d, method="ward.D2")
load_package("sparcl")
png(concat(IMAGES_DIR,'/faculty hierarchical clustering.png'), 
    width = 1024, height = 512)
ColorDendrogram(fit, 
                y = fit5$cluster, 
                main = "Hierarchical Clustering of Faculty Data", 
                xlab = "Euclidean Distance",
                sub = "with Ward D2 Clustering",
                branchlength = 50)
# draw red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")
dev.off()

############################
# answer to c.             #
############################
clvecd <- as.integer(data[,"facrank"])
png(concat(IMAGES_DIR,'/faculty cluster plot3.png'), 
    width = 1024, height = 1024)
plotcluster(x=scaled_data, 
            clvecd=fit5$cluster,
            method="dc",
            clnum=clvecd,
            main="Cluster Plot of Faculty Data with k=5")
dev.off()
png(concat(IMAGES_DIR,'/faculty cluster plot4.png'), 
    width = 1024, height = 1024)
plotcluster(x=scaled_data, 
            clvecd=clvecd,
            method="dc",
            clnum=fit5$cluster,
            main="Cluster Plot of Faculty Data with k=5")
legend("topleft", legend = paste(lda$lev), pch=16, col=1:5)
dev.off()




############################
# example from class       #
############################
load_package("MASS")
set.seed(123456789)
A <- data.frame(source="1", mvrnorm(500, c(0,0), matrix(c(5,0,0,10),2,2)))
B <- data.frame(source="2", mvrnorm(500, c(3,0), matrix(c(3,0,0,1),2,2)))
C <- data.frame(source="3", mvrnorm(500, c(10,10), matrix(c(5,2,2,5),2,2)))
D <- rbind(A,B,C)
plot(D[,2:3], col=D[,1], pch=16)
fit3 <- kmeans(D[,2:3], 3)
plot(D[,2:3], col=fit3$cluster, pch=16)
A_ctr <- fit3$centers[1,]
B_ctr <- fit3$centers[2,]
C_ctr <- fit3$centers[3,]
points(A_ctr["X1"], A_ctr["X2"], pch=4, col="brown", cex=3)
points(B_ctr["X1"], B_ctr["X2"], pch=4, col="brown", cex=3)
points(C_ctr["X1"], C_ctr["X2"], pch=4, col="brown", cex=3)

# hieararchical clustering on sample data
d <- dist(D[,2:3], method="euclidean")
# Option "ward.D2" implements Ward's (1963) clustering criterion 
# (Murtagh and Legendre 2014). With the latter, the dissimilarities 
# are squared before cluster updating.
fit <- hclust(d, method="ward.D2")
load_package("sparcl")
clvecd <- fit3$cluster
ColorDendrogram(fit, 
                y = clvecd,
                main = "Hierarchical Clustering of Multivariate Random Normal Data", 
                xlab = "Euclidean Distance",
                sub = "with Ward D2 Clustering",
                branchlength = 50)
# draw red borders around the 3 clusters 
rect.hclust(fit, k=3, border="red")

# Perform LDA
lda <- lda(formula=D$source ~ ., data=D)
lda$scaling
# Box's M: It performs the Box's M-test for homogeneity of covariance matrices 
# obtained from multivariate normal data according to one classification factor. 
# The test is based on the chi-square approximation.
load_package("biotools")
boxM(D[,2:3], grouping=D[,1])
# Box's M test suggests that the variance-covariance matrices are not
# homegeneous:
# Box's M-test for Homogeneity of Covariance Matrices
# 
# data:  data[, 2:13]
# Chi-Sq (approx.) = 572.96, df = 234, p-value < 2.2e-16

# make predictions
plda <- predict(lda)

# Confusion Matrix:
(confusion_matrix <- table(plda$class, D$source))
# estimate the percentage of faculty rankings that will be mis-classified
round(1 - diag(prop.table(confusion_matrix)), 4)
# total percent incorrect
round(1 - sum(diag(prop.table(confusion_matrix))), 4)
# cross-validation with leave-one-out
lda2 <- lda(D$source ~ ., data=D[,2:3], CV=TRUE)
# Confusion Matrix:
(confusion_matrix2 <- table(lda2$class, D$source))
# estimate the percentage of beetles that will be mis-classified
round(1 - diag(prop.table(confusion_matrix2)), 4)
# total percent incorrect
round(1 - sum(diag(prop.table(confusion_matrix2))), 4)

# plot the LDA projection
(prop.lda = lda$svd^2/sum(lda$svd^2))
lda_data <- data.frame(source = D[,"source"], lda = plda$x)
load_package("ggplot2")
load_package("scales")
ggplot(lda_data) + 
    geom_point(aes(lda.LD1, lda.LD2, col = source, shape = source), size = 2.5) + 
    labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
         y = paste("LD2 (", percent(prop.lda[2]), ")", sep="")) +
    ggtitle("LDA Projection of Random Data") + 
    theme(plot.title = element_text(lineheight=.8, face="bold"))

# Determine number of clusters
wss <- (nrow(D)-1)*sum(apply(D,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(D, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
# K-Means Cluster Analysis for k=3 on original data
fit3 <- kmeans(D, 3) # 3 cluster solution

# cluster plot with ellipses
load_package("cluster")
clusplot(D, 
         fit3$cluster, 
         color=TRUE, 
         shade=TRUE, 
         labels=4, 
         lines=0,
         main="Cluster Plot of Random Data with k=3")

# Centroid Plot against 1st 2 discriminant functions
load_package("fpc")
plotcluster(D[,2:3], 
            method="dc",
            fit3$cluster,
            main="Cluster Plot of Random Data with k=3")
clvecd <- as.integer(D[,1])
plotcluster(x=D, 
            clvecd=fit3$cluster,
            method="dc",
            clnum=clvecd,
            main="Cluster Plot of Random Data with k=3")

plotcluster(x=D[,2:3], 
            clvecd=clvecd,
            method="dc",
            clnum=fit3$cluster,
            main="Cluster Plot of Random Data with k=3")