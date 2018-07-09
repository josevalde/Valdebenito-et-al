# ------------------------------------- 
#loading required packages
# -------------------------------------

library(phytools)
library(ape)
library(metafor)
library(geiger)


# -------------------------------------
# Meta analysis for BLOOD parasites
# -------------------------------------

# -------------------------------------
# loading tree
# -------------------------------------

B.tree <- read.tree("tree blood")

is.binary.tree(B.tree)

B.tree <- multi2di(B.tree)


# -------------------------------------
# loading dataset
# -------------------------------------

B.df <- read.csv("data blood.csv", sep = ";", dec = ",")


# -------------------------------------
# running the model
# -------------------------------------

B.esc <- escalc(ai = M_pos, bi= M_neg, ci=F_pos, di=F_neg, n1i = M_Ex, n2i = F_Ex, 
                data = B.df, measure = "OR") # getting the effect sizes

B.vcvma <- vcv(B.tree, model = "Brownian") # getting variance co-variance matrix

B.m1R <- rma.mv(yi,vi, data = B.esc, method = "ML",
                random = list(~1|species, ~1|method), 
                R=list(species=B.vcvma),
                slab=paste(species, sep = " "))
summary(B.m1R)

ranktest(B.m1R) # checking for publication bias using Kendall's test





# -------------------------------------
# Meta analysis for GASTROINTESTINAL parasites
# -------------------------------------

# -------------------------------------
# loading tree
# -------------------------------------

G.tree <- read.tree("tree GI") #loading tree


# -------------------------------------
# loading dataset
# -------------------------------------

G.df <- read.csv("data GI.csv", sep=";", dec=",")


# -------------------------------------
# running the model
# -------------------------------------

G.esc <- escalc(ai = M_inf, bi= M_not_inf, n1i = M_ex, ci=F_inf, di=F_not_inf, n2i = F_ex, 
                data = G.df, measure = "OR") # getting the effect sizes

G.vcvma <- vcv(G.tree, model = "Brownian") # getting variance co-variance matrix

G.m1 <- rma.mv(yi,vi, data = G.esc, method = "ML",
               random = ~1|species, 
               R=list(species=G.vcvma), 
               slab=paste(species, sep = " ")) #running meta analysis with phylogeny | The slab=... makes the name of the species appear in the final plot
summary(G.m1)

ranktest(G.m1) # checking for publication bias using Kendall's test
