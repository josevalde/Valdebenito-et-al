
# ------------------------------------- 
#loading required packages
# -------------------------------------

library(ape)
library(geiger)
library(caper)


# -------------------------------------
# loading tree
# -------------------------------------

t <- read.nexus("Valdebenito et al. tree")

is.binary.tree(t) # check if tree has polytomies

t <- multi2di(t) # solving polytomies by randomly adding brach lenght of 1e-08 (default). Because is random, you might get differnt topogies 


# -------------------------------------
# loading dataset
# -------------------------------------

data.raw <- read.csv("Valdebenito et al. dataset.csv", sep = ";", dec = ",", header = T) #in this case for name.check you have to put row.names = 2

name.check(t, data.raw) # all good


# -------------------------------------
# normalization of variables
# -------------------------------------

# prevalence blood parasites
data.raw$Male_prev_b <- asin(sqrt((data.raw$Male_prev_b/100))) # males
data.raw$Female_prev_b <- asin(sqrt((data.raw$Female_prev_b/100))) # females

# prevalence gastrointestinal parasites
data.raw$Male_prev_GI <- asin(sqrt((data.raw$Male_prev_GI/100))) # males
data.raw$Female_prev_GI <- asin(sqrt((data.raw$Female_prev_GI/100))) # females

# mortality
data.raw$Male_mort <- asin(sqrt(data.raw$Male_mort)) # males
data.raw$Female_mort <- asin(sqrt(data.raw$Female_mort)) # females

#body mass
data.raw$Male_mass <- log(data.raw$Male_mass) # males
data.raw$Female_mass <- log(data.raw$Female_mass) # females


# -------------------------------------
# creating datasets separate for blood and gastrointestinal parasites
# -------------------------------------

d.b <- data.raw # blood parasites
d.b$Male_prev_GI <- NULL
d.b$Female_prev_GI <- NULL
d.b$Bias_prev_GI <- NULL
View(d.b)
d.g <- data.raw # gastrointestinal parasites
d.g$Male_prev_b <- NULL
d.g$Female_prev_b <- NULL
d.g$Bias_prev_b <- NULL


# -------------------------------------
# Phylogenetic comparative analyses (PGLS)
# -------------------------------------
# Prevalence blood parasites
# -------------------------------------

# blood parasite prevalence ~ body mass. Note that we make a different dataset for every analysis becuause NAs in columns that I am not analyzing deletes the whole row.

d.b1<- d.b[,c(2,9,6)] # PGLS MALES
comp1 <- comparative.data(t, d.b1, vcv=T , species,vcv.dim = 3)
m1 <- pgls(Male_prev_b ~ Male_mass, comp1, lambda = "ML")
summary(m1) 
plot(Male_prev_b ~ Male_mass, data=d.b1)
abline(m1)

d.b2<- d.b[,c(2,7,10)] # PGLS FEMALES
comp2 <- comparative.data(t, d.b2, vcv=T , species,vcv.dim = 3)
m2 <- pgls(Female_prev_b ~ Female_mass, comp2, lambda = "ML") 
summary(m2)
plot(Female_prev_b ~ Female_mass, data=d.b2)
abline(m2)

d.b3 <- d.b[,c(2,8,11)] #PGLS BIAS
comp3 <- comparative.data(t, d.b3, vcv=T , species,vcv.dim = 3)
m3 <- pgls(Bias_prev_b ~ SSD, comp3, lambda = "ML")
summary(m3)
plot(Bias_prev_b ~ SSD, data = d.b3)
abline(m3)


# -------------------------------------
# Prevalence gastrointestinal parasites
# -------------------------------------

# gastrointestinal parasite prevalence ~ body mass

d.g1<- d.g[,c(2,9,6)] # PGLS MALES
comp4 <- comparative.data(t, d.g1, vcv=T , species,vcv.dim = 3)
m4 <- pgls(Male_prev_GI ~ Male_mass, comp4, lambda = "ML")
summary(m4) 
plot(Male_prev_GI ~ Male_mass, data=d.g1)
abline(m4)

d.g2<- d.g[,c(2,10,7)] # PGLS FEMALES
comp5 <- comparative.data(t, d.g2 , vcv=T, species,vcv.dim = 3)
m5 <- pgls(Female_prev_GI ~ Female_mass, comp5, lambda = "ML") 
summary(m5)
plot(Female_prev_GI ~ Female_mass, data=d.g2)
abline(m5)

d.g3 <- d.g[,c(2,8,11)] # PGLS BIAS
comp6 <- comparative.data(t, d.g3, vcv=T , species,vcv.dim = 3)
m6 <- pgls(Bias_prev_GI ~ SSD, comp6, lambda = "ML")
summary(m6)
plot(Bias_prev ~ SSD, data = d.g3)
abline(m6)


# -------------------------------------
# mortalities & blood parasites
# -------------------------------------

#mort ~ prev blood + body mas + poly score

d.b4<- d.b[,c(2,3,6,9,12)] # PGLS MALES
comp7 <- comparative.data(t, d.b4, species, vcv = T, vcv.dim = 3)
m7 <- pgls(Male_mort ~ Male_prev_b + Male_mass + Male_poly, comp7, lambda = "ML")
summary(m7)
summary(d.b4)
d.b5<- d.b[,c(2,4,7,10,13)] # PGLS FEMALES
comp8 <- comparative.data(t, d.b5, species, vcv = T, vcv.dim = 3)
m8 <- pgls(Female_mort ~ Female_prev_b + Female_mass + Female_poly, comp8, lambda = "ML")
summary(m8)

d.b6<- d.b[,c(2,5,8,11,14)] # PGLS BIAS
comp9 <- comparative.data(t, d.b6, species, vcv = T, vcv.dim = 3)
m9 <- pgls(Bias_mort ~ Bias_prev_b + SSD + Bias_poly, comp9, lambda = "ML")
summary(m9)

# mort ~ prev blood

d.b7<- d.b[,c(2,3,9)] # PGLS MALES
comp10 <- comparative.data(t, d.b7, species, vcv = T, vcv.dim = 3)
m10 <- pgls(Male_mort ~ Male_prev_b, comp10, lambda = "ML")
summary(m10)

d.b8<- d.b[,c(2,4,10)] # PGLS FEMALES
comp11 <- comparative.data(t, d.b8, species, vcv = T, vcv.dim = 3)
m11 <- pgls(Female_mort ~ Female_prev_b, comp11, lambda = "ML")
summary(m11)

d.b9<- d.b[,c(2,5,11)] # PGLS BIAS
comp12 <- comparative.data(t, d.b9, species, vcv = T, vcv.dim = 3)
m12 <- pgls(Bias_mort ~ Bias_prev_b, comp12, lambda = "ML")
summary(m12)


# -------------------------------------
# mortalities & gastrointestinal parasites
# -------------------------------------

#mort ~ prev GI + body mas + poly score

d.g4<- d.g[,c(2,3,6,9,12)] # PGLS MALES
comp13 <- comparative.data(t, d.g4, species, vcv = T, vcv.dim = 3)
m13 <- pgls(Male_mort ~ Male_prev_GI + Male_mass + Male_poly, comp13, lambda = "ML")
summary(m13)

d.g5<- d.g[,c(2,4,7,10,13)] # PGLS FEMALES
comp14 <- comparative.data(t, d.g5, species, vcv = T, vcv.dim = 3)
m14 <- pgls(Female_mort ~ Female_prev_GI + Female_mass + Female_poly, comp14, lambda = "ML")
summary(m14)

d.g6<- d.g[,c(2,5,8,11,14)] # PGLS BIAS
comp15 <- comparative.data(t, d.g6, species, vcv = T, vcv.dim = 3)
m15 <- pgls(Bias_mort ~ Bias_prev_GI + SSD + Bias_poly, comp15, lambda = "ML")
summary(m15)


#mort ~ prev GI

d.g7<- d.g[,c(2,3,9)] # PGLS MALES
comp16 <- comparative.data(t, d.g7, species, vcv = T, vcv.dim = 3)
m16 <- pgls(Male_mort ~ Male_prev_GI, comp16, lambda = "ML")
summary(m16)

d.g8 <- d.g[,c(2,4,10)] # PGLS FEMALES
comp17 <- comparative.data(t, d.g8, species, vcv = T, vcv.dim = 3)
m17 <- pgls(Female_mort ~ Female_prev_GI, comp17, lambda = "ML")
summary(m17)

d.g9<- d.g[,c(2,5,11)] # PGLS BIAS
comp18 <- comparative.data(t, d.g9, species, vcv = T, vcv.dim = 3)
m18 <- pgls(Bias_mort ~ Bias_prev_GI, comp18, lambda = "ML")
summary(m18)
