
library(dplyr)
library(ggplot2)
library(scatterpie)
library(grafify)
library(foreach)
library(doParallel)
library(cluster)
library(ape)
library(BAT)
library(cowplot)
library(hypervolume)
library(VGAM)
library(rlang)
library(data.table)
library(picante)
library(gridExtra)
library(grid)
library(tidyr) #needs to >= version 1.3.1
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")

##########################################################################
### Load and Format Data############################################
##############################################################################

#Load in AVONET Trait data
extant <- read.csv("Data\\AVONET_BIRDLIFE.csv")

#load in IUCN and island endemic classifications
sayol <- read.csv("Data\\allSp_birdlife.csv")

#if seabird sensitivity test, swap island endemic columns
if (seabird){
  sayol$IslandEndemic <- sayol$IslandEndemic2
}

#Load in extinct species data
extinct <- read.csv("Data\\ExtinctImputedTraits_Averages.csv")

#Change ExtinctionType and ExtinctionPeriod codes (these were updated in the
#datasets during review).
extinct$ExtinctionType <- sapply(extinct$ExtinctionType,
                                 function(y) {switch(y,
                                                     "EX_U" = "EPU",
                                                     "EX_A" = "EPA",
                                                     "EX_IUCN" = "IUCN",
                                                     y)})

extinct$ExtinctionPeriod <- sapply(extinct$ExtinctionPeriod,
                                   function(y) {switch(y,
                                                       "EX_Other" = "All",
                                                       "EX_IUCN" = "IUCN",
                                                       y)})

#create HWI column
extinct <- mutate(extinct, Hand.Wing.Index = 100 * (Kipps.Distance / Wing.Length))

if (length(sayol$species[which(!sayol$species %in% extant$species)]) != 0) stop("wishbone") 

if (length(extant$species[which(!extant$species %in% sayol$species)]) != 6) stop("wishbone A") 

##check species in both and remove:
#there are six species listed in both. Five are EW and we will class these as extinct
#The 6th (Zosterops conspicillatus) has recently been split and the nominate subpsecies
#i.e. this name, is now considered extinct with the other subspecies extant and given a new
#name (Z. saypani). As such, remove all six from AVONET.
wea3 <- which(extant$species %in% extinct$species)

extant <- extant[-wea3,]
if (!identical(sort(extant$species), sort(sayol$species))) stop("wishbone B")

#change EW to IUCN
if (!identical(which(extinct$ExtinctionPeriod == "EW"), 
               which(extinct$ExtinctionType == "EW"))){
  stop("moonshadow")
}
extinct$ExtinctionPeriod[which(extinct$ExtinctionPeriod == "EW")] <- "IUCN"
extinct$ExtinctionType[which(extinct$ExtinctionType == "EW")] <- "IUCN"

#add in extant status and ExtinctionPeriod Columns
extant$Extant <- "Extant"
extant$ExtinctionPeriod <- "Current"

#take island endemic column from sayol and add into AVONET
extant2 <- full_join(extant, sayol, by = "species")
if(!nrow(extant) == nrow(extant2)) stop("Error 1")
extant2 <- dplyr::select(extant2, -c(Order.y, Family.y, Genus.y))

extinct$Extant <- "Extinct"

#ExtinctionType is the column which distinguishes between unknown (EPU)
#and anthropogenic (EPA) extinctions. For current species, just say "Current"
extant2$ExtinctionType <- "Current"

##select required columns
extant3 <- dplyr::select(extant2, species, "Genus" = Genus.x, 
                         "Family" = Family.x, "Order" = Order.x, Extant, 
                         ExtinctionPeriod, ExtinctionType,
                         IslandEndemic, status, Hand.Wing.Index,
                         Beak.Length_Culmen,	Beak.Length_Nares,	Beak.Width,
                         Beak.Depth,	Tarsus.Length,	Wing.Length,
                         Kipps.Distance,	Tail.Length,	Mass)

extinct$status <- "EX"

extinct2 <- dplyr::select(extinct, species, Genus, 
                          Family, Order, Extant, ExtinctionPeriod, 
                          ExtinctionType, IslandEndemic,
                          status, Hand.Wing.Index,
                          Beak.Length_Culmen,	Beak.Length_Nares,	Beak.Width,
                          Beak.Depth,	Tarsus.Length,	Wing.Length,
                          Kipps.Distance,	Tail.Length,	Mass)

#if doing the individual imputation sensitivity test,
#replace the 9 morpho traits in extinct2 with the values
#from an individual run (ind_impu)
if (!is.null(ind_impu)){
  
  if (!all(extinct2$species %in% ind_impu$species)){
    stop("mardy")
  }
  
  ind_impu_CN <- c("Beak.Length_Culmen", "Beak.Length_Nares",
                   "Beak.Width", "Beak.Depth", "Tarsus.Length",
                   "Wing.Length", "Kipps.Distance", "Tail.Length",
                   "Mass")
  
  if (!all(ind_impu_CN %in% colnames(datCom[[i]]))){
    stop("monkey")
  }
  
  ind_impu2 <- filter(ind_impu, species %in% extinct2$species)
  
  #match row orders
  ind_impu_M <- match(extinct2$species, ind_impu2$species)
  ind_impu2 <- ind_impu2[ind_impu_M,]
  if (!identical(ind_impu2$species, extinct2$species)){
    stop("rosefinch")
  }
  if (!nrow(ind_impu2) == 610) stop("arctic")
  
  extinct2[,ind_impu_CN] <- ind_impu2[,ind_impu_CN]
  
  ppp = paste0(ceiling(runif(1,1,1000)),".csv")
  write.csv(extinct2$Beak.Width, file = ppp)
  
}#eo ind_impu


##merge datasets

if(!identical(colnames(extinct2), colnames(extant3))) stop("hot dog")

allSp <- bind_rows(extant3, extinct2)


####load in future simulated extinction data
future <- readRDS("Data\\future.rds")

####Phylogeny#######################################################

if (FD_PD == "PD"){

if (seabird){
  future <- future[which(names(future) %in% c("Birdtree_all", "Birdtree_isl2"))]
} else {
  future <- future[which(names(future) %in% c("Birdtree_all", "Birdtree_isl"))]
}

phylo_all <- readRDS("Data\\JetzExtinct50Trees.rds")

#this is only for checks, as extinct2/3 not used again for phylo analyses
extinct3 <- extinct2

if (!all(extinct3$species %in% phylo_all[[1]]$tip.label)){
  stop("Pneuma")
}

#Note there are 12 species in our extinct dataset that are in Jetz
if(!length(phylo_all[[1]]$tip.label) == 9993 + nrow(extinct2) - 12){
  stop("scar tissue")
}

##load allSp version for phylo analyses
allSp_phylo <- read.csv("Data\\allSp_phylo.csv")

#Change ExtinctionType and ExtinctionPeriod codes (these were updated in the
#datasets during review)
allSp_phylo$ExtinctionType <- sapply(allSp_phylo$ExtinctionType,
                                     function(y) {switch(y,
                                                         "EX_U" = "EPU",
                                                         "EX_A" = "EPA",
                                                         "EX_IUCN" = "IUCN",
                                                         y)})

allSp_phylo$ExtinctionPeriod <- sapply(allSp_phylo$ExtinctionPeriod,
                                       function(y) {switch(y,
                                                           "EX_Other" = "All",
                                                           "EX_IUCN" = "IUCN",
                                                           y)})


#if seabird sensitivity test, swap island endemic columns
if (seabird){
  allSp_phylo$IslandEndemic <- allSp_phylo$IslandEndemic2
}

colnames(allSp_phylo)[which(colnames(allSp_phylo) == "IUCN")] <-
  "status"

if (!identical(sort(filter(allSp_phylo, Extant == "Extinct")$species),
               sort(extinct3$species)) | 
    !identical(sort(filter(allSp_phylo, Extant == "Extinct")$species),
               sort(extinct3$species))){
  stop("moonshadow")
}

if (!identical(nrow(allSp_phylo), length(phylo_all[[1]]$tip.label))){
  stop("Father and son")
}

if(!all(allSp_phylo$species %in% phylo_all[[1]]$tip.label)){
  stop("the wind")
}

if (!identical(sort(extinct3$species[which(extinct3$ExtinctionType == "EPU")]), 
               sort(allSp_phylo$species[which(allSp_phylo$ExtinctionType == "EPU")]))){
  stop("moonshadow2")
}

allSp_phylo$ExtinctionPeriod[which(allSp_phylo$ExtinctionPeriod == "EW")] <- "IUCN"
allSp_phylo$ExtinctionType[which(allSp_phylo$ExtinctionType == "EW")] <- "IUCN"

allSp2 <- allSp_phylo
rownames(allSp2) <- allSp_phylo$species

######################################################################
## format traits;

} else if (FD_PD == "FD"){
  
if (seabird){
  future <- future[which(names(future) %in% c("Birdlife_all", "Birdlife_isl2"))]
} else {
  future <- future[which(names(future) %in% c("Birdlife_all", "Birdlife_isl"))]
}

  traits <-
    dplyr::select(
      allSp,
      "Beak.Length_Culmen",
      "Beak.Length_Nares",
      "Beak.Width",
      "Beak.Depth",
      "Tarsus.Length",
      "Wing.Length",
      "Kipps.Distance",
      "Tail.Length",
      "Mass"
    )
  rownames(traits) <- allSp$species

##KIWIs & Kipps################################
#the kiwi species are extreme outliers in regard to wing length and
#tail length, but not Kipps (as other species have Kipps of 0.1). 
#What we will do is replace them with the mean non-kiwi 
#observed trait values (across extant and extinct).
w3 <- which(rownames(traits) %in% c("Apteryx_australis", "Apteryx_haastii",
                                    "Apteryx_mantelli","Apteryx_owenii",
                                    "Apteryx_rowi"))
#create version without kiwis
traits_noKiwi <- traits[-w3,]

#work out smallest trait values for remaining species (or mean)
mean_noKiwi <- apply(traits_noKiwi[,c("Wing.Length",
                                      "Tail.Length")], 2, mean)

#switch kiwi values for these min values
cn_noKiwi <- c("Wing.Length", "Tail.Length")
for (i in 1:2){
  traits[w3,cn_noKiwi[i]] <- mean_noKiwi[i]
}
#check
if (!identical(apply(traits[,c("Wing.Length", "Tail.Length")], 2, mean),
               mean_noKiwi)) stop("Township rebellion")

#Now need to do the same for the moa
wM <- filter(extinct, Order == "Dinornithiformes")$species

wM3 <- which(rownames(traits) %in% wM)
if (!all(traits[wM3,"Wing.Length"] == 0.1)){
  stop("Moa")
}
traits[wM3,"Wing.Length"] <- mean_noKiwi[1]
#HWI is in allSp
if (!all(allSp[wM3,"Hand.Wing.Index"] == 100)){
  stop("Moa2")
}
allSp[wM3, "Hand.Wing.Index"] <- 0.1

##For the non-measurable Kipp's and HWI species, they
#all have Kipp's and HWI < 1. 
#To reduce the outlier effect, replace with the value 1.
#Kipps inside the traits df; also in allSp but those columns
#not used downstream, so no need to change.
traits_kipps <- which(traits$Kipps.Distance < 1)
#rownames(traits)[traits_kipps]#should be 29
if (length(traits_kipps) != 29 | 
    !all(traits$Kipps.Distance[traits_kipps] < 1)) stop("Parker")
traits$Kipps.Distance[traits_kipps] <- 1

#HWI is not stored in traits, but in original allSp df
allSp_hWI <- which(allSp$Hand.Wing.Index < 1)
#allSp$species[allSp_hWI]#should be 29
if (length(allSp_hWI) != 29 | 
    !all(allSp$Hand.Wing.Index[allSp_hWI] < 1)) stop("Hanson")
allSp$Hand.Wing.Index[allSp_hWI] <- 1


traits1 <- apply(traits, 2, log) %>% as.data.frame() #log-transform following Pigot et al

###########################################################
##USE BODY-SIZE CORRECTED TRAITS OR NOT
#body_size_CORR <- FALSE
#body_size_CORR <- TRUE #use BS corrected traits
###############################################################

if (body_size_CORR){
  
  ######body mass regression approach 
  #function to return residuals from trait-mass regression
  #resp = morpho trait to use as a response
  trait_residuals <- function(traitExtant, resp){
    modBlur <- lm(traitExtant[,resp] ~ traitExtant$Mass)
    residBlur <- residuals(modBlur)
    return(residBlur)
  }
  
  residTraits <- vapply(c("Beak.Length_Culmen", "Beak.Length_Nares", 
                          "Beak.Width", 
                          "Beak.Depth", "Tarsus.Length",
                          "Wing.Length", "Kipps.Distance","Tail.Length"), 
                        function(x) trait_residuals(traits1, resp = x),
                        FUN.VALUE = numeric(nrow(traits1)))
  
  residTraits <- cbind(residTraits, "Mass" = traits1$Mass) 
  
  traits2 <- apply(residTraits, 2, scale)#we are scaling traits prior to PCA
  rownames(traits2) <- rownames(traits)
  
} else {
  
  traits2 <- apply(traits1, 2, scale)#we are scaling traits prior to PCA
  rownames(traits2) <- rownames(traits)

}

#if running the sensitivity test removing individual traits
if (!is.null(tra_remove)){
  traits2 <- traits2[,-which(colnames(traits2) == tra_remove)]
  if (ncol(traits2) != 8) stop("the next episode")
}

trai.pca <- prcomp(traits2)
traits3 <- trai.pca$x #all nine PCA axes

#If running the body shape analysis, we need to remove
#the first PCA axis
if (body_shape){
  traits3_CN <- colnames(traits3)
  traits3 <-  traits3[,2:ncol(traits3)]
  colnames(traits3) <- traits3_CN[1:8]
}


#if doing the trait removal sensitivty test, can't do the beak analyses
#as one of the beak traits may be removed
if (is.null(tra_remove)){

#do beak specific PCA
PCA_beak <- prcomp(traits2[,c("Beak.Length_Culmen",	"Beak.Length_Nares",
                              "Beak.Width",	"Beak.Depth")])#first 3 axis = 99% of variation
PCA_beaks <- PCA_beak$x
colnames(PCA_beaks) <- c("PC1_beak", "PC2_beak", "PC3_beak", "PC4_beak")

######dendrogram for all species using nj method in BAT
load("Data\\dendro_all.RData")
if (!identical(sort(rownames(traits3)), sort(dendro_all$tip.label))) stop("Dylan")
if (!identical(sort(allSp$species), sort(dendro_all$tip.label))) stop("Dylan")

# create main allSp2 dataframe
allSp2 <-  cbind(allSp, traits3, PCA_beaks) %>% as.data.frame()

} else {
  load("Data\\dendro_all.RData")
  if (!identical(sort(rownames(traits3)), sort(dendro_all$tip.label))) stop("Dylan")
  if (!identical(sort(allSp$species), sort(dendro_all$tip.label))) stop("Dylan")
  allSp2 <-  cbind(allSp, traits3) %>% as.data.frame()
}#eo if tra_remove


################################################################
############Hypervolume Default parameters#####################
################################################################

#These are the functions inside the hypervolume package, which
#calculate the argument based on the data provided. Here, we calculate
#them for the full dataset (allSp2) to ensure consistent values
#used across analyses.

##samples.per.point
#we use 5 PCA axes in the main analyses;
#the beak analyses only use 3, but we still use this value to be consistent
SPP <- ceiling((10^(3 + sqrt(no_axes)))/nrow(allSp2))

##bandwidth - only for GAUSSIAN
#BAND_p needed to extract the correct PCA axes from allSp2
BAND_p <- paste0("PC", no_axes)
BAND <- hypervolume::estimate_bandwidth(dplyr::select(allSp2, 
                                                      "PC1":all_of(BAND_p)))
###################################################################

} else{
  stop ("FD_PD should be one of FD or PD")
}

prehistoric_all <- allSp2 %>% mutate("Period" = "Prehistoric:all")

#######################################################################
##########FUNCTIONS###########################################
#############################################################

##Percentage change in Alpha diversity between different time period:
#P2C, H2C, C2F
perc_change <- function(y){
  
  perc_change_int <- function(v1,v2) ((v2 - v1) / v1) * 100  #v1 = original value
  
  r1 <- perc_change_int(v1 = y[1], v2 = y[2])#P2h 
  r2 <- perc_change_int(v1 = y[1], v2 = y[3])#P2c
  r3 <- perc_change_int(v1 = y[2], v2 = y[3])#H2C
  r4 <- perc_change_int(v1 = y[3], v2 = y[4])
  
  return(c(r1,r2,r3,r4))
}

#dat should be the prehistoric or historic dataset versions,
#n = number of extinctions
#trait = trait to focus on
ind_null <- function(dat, trait, n){
  
  #random sample of n sp to go extinct
  n2 <- sample(1:nrow(dat), n, replace = F)
  #remove extinct sp
  dat2 <- dat[-n2,]
  #calculate median and sd of trait
  null_trait <- dat2[,trait]
  r1 <- median(null_trait)
  r2 <- sd(null_trait)
  r3 <- c(r1, r2)
  return(r3)
}


##Get Z-scores for each null distribution (and 2-tailed P-value)
#if null distributions roughly normal use standard P-value and SES,
#if not use the Lhotsky ES (type == "ES") and 2-sided p-value 

#NOTE that for ES, the maximum value you can seemingly get is an
#ES value of 3.29, when you have 999 null iterations. When n = 999,
#if all the null values are less than the obs, it returns a probability
#of 0.9995 which the probit of is 3.29. Increasing n to 9999 increase the 
#prob to 0.99995 and ES to 3.89, and so on. 

#for the ES values the p-value is on a continuous scale,
#from 0-1, with values <0.025 and >0.975 being significant at the 0.05
#level. Whereas for the SES, the probability is for the lower tail which
#is then multiplied by 2 to give you the two-tailed test, where 
#any P < 0.05 is significant. It needs to be divided by 2 if you want
#to directly compare with the ES version

zP <- function(dis, obs, type = "ES"){
  if (type == "SES"){
    z <- (obs - mean(dis)) / sd(dis)
    p <- 2*pnorm(-abs(z))
  } else if (type =="ES"){
    dis <- c(obs, dis)# don't forget to add the obs into the null values
    res <- (sum(dis < obs) + sum(dis == obs)/2) / length(dis)
    z <- VGAM::probitlink(res)
    p <- res
  }
  return(c(z,p))
}

####Individual trait function: mass and Kipp's


#prehistoric_all = all species across all time periods (should = prehistoric_all)
#extinct = main extinct dataset
#future_F = one run of future extinctions from the
#main future dataset (from iucn_sim)
##Set TRAIT to one of Mass or Hand.Wing.Index
##set GEOG to one of All, "Yes (for Isl endemic)
##Set TEST to one of Median or SD
#n = number of null iterations
#titl = ggplot title
#bandwidth = geom violin bw argument
#color = violin colors
#hist = if true plots histograms of null values instead (3 overlaid)
#analysis = main or beak
#beak_data = data for doing the beak analysis

null_plot <- function(prehistoric_all_F = prehistoric_all,
                      extinct_F = extinct,
                      future_F = future,
                      TRAIT, GEOG, TEST, n = 99, titl = "a)", 
                      bandwidth = 9999,
                      #  color = "blue", 
                      hist = FALSE, 
                      analysis = "main", beak_data = NULL){
  
  if (analysis == "main"){
    #this then filters or sets the dataset to the correct version based on selected
    #GEOG value
    if (GEOG == "All"){
      dat_ind <- prehistoric_all_F
    } else {
      dat_ind <- filter(prehistoric_all_F, IslandEndemic == GEOG)
    }
    
    ##extinction numbers
    Hist_N <- table(dat_ind$ExtinctionPeriod)["IUCN"] %>% as.vector()
    Pre_N <- table(dat_ind$ExtinctionPeriod)["All"] %>% as.vector()
    All_N <- Hist_N + Pre_N
    if (GEOG == "All"){
      if (All_N != nrow(extinct_F)) stop ("error 5")
    }
    Fut1_N <- length(future_F[[1]])#all elements have same N
    
    #get obs median trait value for prehistoric, historit, current and 2 future periods
    obs_P <- dat_ind %>%
      dplyr::select(all_of(TRAIT)) %>%
      summarize(Median = median(!!rlang::sym(TRAIT)), #uses rlang to allow a character be used instead of colname
                SD = sd(!!rlang::sym(TRAIT)))
    
    obs_H <- filter(dat_ind, ExtinctionPeriod %in% c("IUCN", "Current")) %>%
      dplyr::select(all_of(TRAIT)) %>%
      summarize(Median = median(!!rlang::sym(TRAIT)), #uses rlang to allow a character be used instead of colname
                SD = sd(!!rlang::sym(TRAIT)))
    
    obs_C <- filter(dat_ind, ExtinctionPeriod == "Current") %>%
      dplyr::select(all_of(TRAIT)) %>%
      summarize(Median = median(!!rlang::sym(TRAIT)), #uses rlang to allow a character be used instead of colname
                SD = sd(!!rlang::sym(TRAIT)))
    
    if (!all(unlist(future_F) %in% 
             filter(dat_ind, ExtinctionPeriod == "Current")$species)){
      stop("cigs")
    }
    
    ##run it across each of the iucn_sim runs; here we want to filter out
    #the species listed in future_F (as these are the extinct sp.)
    obs_F1 <- sapply(future_F, function(x){
      f1b <- filter(dat_ind, ExtinctionPeriod == "Current",
                    !species %in% x) %>%
        dplyr::select(all_of(TRAIT)) %>%
        summarize(Median = median(!!rlang::sym(TRAIT)), #uses rlang to allow a character be used instead of colname
                  SD = sd(!!rlang::sym(TRAIT)))
      unlist(f1b)
    })
    
    if (!is.matrix(obs_F1) || 
        ncol(obs_F1) != length(future_F) ||
        !identical(rownames(obs_F1),c("Median", "SD"))){
      stop("alcohol")
    }
    
    obs_F <- apply(obs_F1, 1, median)#take median values
    
    #merge (obs_C used twice as we have two comparisons using current period)
    obs_null <- bind_rows(obs_C, obs_C, obs_F) %>%
      mutate(Type = c("All", "IUCN", "Fut")) 
    colnames(obs_null)[1:2] <- c("Median", "SD")
    obs_null$Type <- factor(obs_null$Type, 
                            levels = c("All", "IUCN", "Fut"))
    
    P2C <- replicate(n, ind_null(dat = dat_ind, trait = TRAIT, n = All_N)) %>%
      t()%>%
      as.data.frame()
    P2C$V3 <- "All"
    
    #for Hist to current comparison, need to get rid of pre-historic species
    datH2c <- filter(dat_ind, ExtinctionPeriod %in% c("IUCN", "Current"))
    H2C <- replicate(n, ind_null(dat = datH2c, trait = TRAIT, n = Hist_N)) %>%
      t()  %>%
      as.data.frame()
    H2C$V3 <- "IUCN"
    
    #for C2F comparison, need to get rid of pre-historic & Historic species
    datC2F1 <- filter(dat_ind, ExtinctionPeriod %in% c("Current"))
    C2F <- replicate(n, ind_null(dat = datC2F1, trait = TRAIT, n = Fut1_N)) %>%
      t()  %>%
      as.data.frame()
    C2F$V3 <- "Fut"
    
    #plot ylabs
    if (TRAIT == "Mass" & TEST == "Median"){
      yl <- "Body Mass (Median)"
    } else if (TRAIT == "Mass" & TEST == "SD"){
      yl <- "Body Mass (SD)"
    } else if (TRAIT == "Log_Mass" & TEST == "Median"){
      yl <- "Log BM (Median)"
    } else if (TRAIT == "Log_Mass" & TEST == "SD"){
      yl <- "Log BM (SD)" 
    } else if (TRAIT == "Kipps.Distance" & TEST == "Median"){
      yl <- "Kipps Distance (Median)"
    } else if (TRAIT == "Kipps.Distance" & TEST == "SD"){
      yl <- "Kipps Distance (SD)"
    } else if (TRAIT == "Log_Kipps.Distance" & TEST == "Median"){
      yl <- "Log KD (Median)"
    } else if (TRAIT == "Log_Kipps.Distance" & TEST == "SD"){
      yl <- "Log KD (SD)" 
    } else if (TRAIT == "Hand.Wing.Index" & TEST == "Median"){
      yl <- "HWI (Median)"
    } else if (TRAIT == "Hand.Wing.Index" & TEST == "SD"){
      yl <- "HWI (SD)"
    } else if (TRAIT == "Log_Hand.Wing.Index" & TEST == "Median"){
      yl <- "Log HWI (Median)"
    } else if (TRAIT == "Log_Hand.Wing.Index" & TEST == "SD"){
      yl <- "Log HWI (SD)"   
    } else{
      yl = "wooooh"
    }
    
    #merge null model runs
    null_all <- bind_rows(P2C, H2C, C2F)
    colnames(null_all) <- c("Median", "SD", "Type")
    null_all$Type <- factor(null_all$Type, 
                            levels = c("All", "IUCN", "Fut"))
    
    
    #set point colour based on significance of ES
    beak_col <- vector(length = length(obs_null$Type))
    
    beak_ZP_mat <- matrix(ncol = 4, nrow = length(obs_null$Type))
    colnames(beak_ZP_mat) <- c("ES", "ES_P", "SES", "SES_P")
    rownames(beak_ZP_mat) <- obs_null$Type
    
    for (i in 1:length(obs_null$Type)){
      
      z <- obs_null$Type[i]
      dum_null_trait <- filter(null_all, Type == z)
      if (nrow( dum_null_trait) != n) stop("error 8")
      null_vals <- dum_null_trait[, TEST]
      obs_val <- obs_null[i, TEST]
      beak_ZP <- zP(null_vals, obs_val, type = "ES")
      beak_ZP_mat[i,1:2] <- beak_ZP
      beak_ZP_mat[i,3:4] <- zP(null_vals, obs_val, type = "SES")
      #Kipps/HWI is two-tailed test, but others one-tailed
      if (TRAIT %in% c("Kipps.Distance", "Log_Kipps.Distance",
                       "Hand.Wing.Index", "Log_Hand.Wing.Index")){
        beak_col[i] <- ifelse(beak_ZP[2] > 0.975 | beak_ZP[2] < 0.025,
                              "#000000", "#ed3832")
      } else{
        beak_col[i] <- ifelse(beak_ZP[2] < 0.05,
                              "#000000", "#ed3832")
      }
    }
    
    beak_ZP_mat <- round(beak_ZP_mat, 2)
    
    #set bandwidth to nrd method if no value given
    if (bandwidth == 9999) bandwidth <- "nrd"
    
    #plot histogram of null values instead if hist == TRUE
    if (hist){
      joris <- ggplot(data = null_all, 
                      aes(!!rlang::sym(TEST))) + 
        theme_bw() + geom_histogram(alpha=0.4, bins = 30) + 
        ggtitle(titl) + facet_wrap(~Type, scales = "free")
    } else {
      
      jorisA <- ggplot(data = null_all, 
                       aes(Type, !!rlang::sym(TEST))) + 
        theme_bw() + xlab("Time Period") + ylab(yl) + 
        geom_violin(aes(fill = Type), fill = "#48ABAE",
                    bw = bandwidth)  +
        scale_fill_grafify(palette = "fishy") +
        geom_point(data = obs_null, aes(Type, !!rlang::sym(TEST)), 
                   size = 5,
                   shape = 18, alpha = 0.6,
                   col =  beak_col)  + ggtitle(titl) + 
        theme(axis.text = element_text(size = 13), 
              axis.title = element_text(size = 15),
              plot.title = element_text(size = 17)) +
        theme(legend.position = "none") + 
        scale_x_discrete(labels=c("All" = expression(All%->%Cur),
                                  "IUCN" = expression(IUCN%->%Cur),
                                  "Fut" = expression(Cur%->%Fut)))
      
      ##make small version just of the three observed values for inset plot
      obs_null2 <- bind_rows(obs_P, obs_H, obs_C, obs_F) %>%
        mutate(Type = c("All", "IUCN", "Cur", "Fut")) 
      colnames(obs_null2)[1:2] <- c("Median", "SD")
      obs_null2$Type <- factor(obs_null2$Type, 
                               levels = c("All", "IUCN", "Cur", "Fut"))
      
      #y-axis lims
      IT_YL <- range(select(obs_null2, !!rlang::sym(TEST)))
      IT_YL_U <- IT_YL[2] + (abs(diff(IT_YL)) * 0.1)
      IT_YL_L <- IT_YL[1] - (abs(diff(IT_YL)) * 0.1)
      
      jorisB <- ggplot(data = obs_null2, 
                       aes(Type, !!rlang::sym(TEST))) + 
        theme_bw() + xlab("") + ylab("")  +
        geom_point(size = 3,
                   shape = 18, col = "#716aa3") + 
        geom_line(data=subset(obs_null2, 
                              Type %in% c("All", "IUCN", "Cur")),
                  group = 1, col = "#716aa3") + 
        geom_line(data=subset(obs_null2, 
                              Type %in% c("Cur","Fut")),
                  group = 1, linetype = 2,
                  col = "#716aa3") +
        ylim(IT_YL_L, IT_YL_U)
      
      joris <- list(list(jorisA, jorisB), obs_null2, beak_ZP_mat, obs_F1)
      
    }
    
  } else if (analysis == "beak"){
    
    obs <- beak_data[[1]] 
    obs$Type <- rownames(beak_data[[1]])
    obs$Type <- factor(obs$Type, levels = c("All", "IUCN", "Cur", "Fut"))
    
    obs_null <- data.frame("Alpha" = c(obs["Cur","Alpha"], 
                                       obs["Cur","Alpha"], obs["Fut","Alpha"]),
                           Type = c("All", "IUCN", "Fut")) 
    obs_null$Type <- factor(obs_null$Type, 
                            levels = c("All", "IUCN", "Fut"))
    
    
    null_all <- data.table::rbindlist(beak_data[[2]], 
                                      fill=FALSE, idcol=NULL)
    
    null_all$Type <- rep(c("All", "IUCN", "Fut"), (nrow(null_all) / 3))
    null_all$Type <- factor(null_all$Type, 
                            levels = c("All", "IUCN", "Fut"))
    
    
    if (hist){
      joris <- ggplot(data = null_all, aes(Alpha, fill = Type)) + 
        theme_bw() + geom_histogram(alpha=0.4, position="identity") + 
        ggtitle(titl) + facet_wrap(~Type, scales = "free")
    } else {
      
      #set point colour based on significance of ES
      beak_col <- vector(length = length(obs_null$Type))
      
      beak_ZP_mat <- matrix(ncol = 4, nrow = length(obs_null$Type))
      colnames(beak_ZP_mat) <- c("ES", "ES_P", "SES", "SES_P")
      rownames(beak_ZP_mat) <- obs_null$Type
      
      for (i in 1:length(obs_null$Type)){
        
        z <- obs_null$Type[i]
        dum_null_trait <- filter(null_all, Type == z)
        null_vals <- dum_null_trait[, "Alpha"]
        obs_val <- obs_null[i, "Alpha"]
        beak_ZP <- zP(null_vals$Alpha, obs_val, type = "ES")
        beak_ZP_mat[i,1:2] <- beak_ZP
        beak_ZP_mat[i,3:4] <- zP(null_vals$Alpha, obs_val, type = "SES")
        beak_col[i] <- ifelse(beak_ZP[2] < 0.05,
                              "#000000", "#ed3832")
      }
      
      beak_ZP_mat <- round(beak_ZP_mat, 2)
      
      jorisA <- ggplot(data = null_all, aes(Type, Alpha)) + 
        theme_bw() + xlab("Time Period") + ylab("Alpha") + 
        geom_violin(aes(fill = Type), fill = "#48ABAE",
                    bw =  "nrd")  +
        scale_fill_grafify(palette = "fishy") +
        geom_point(data = obs_null, 
                   aes(Type, Alpha), 
                   size = 5,
                   shape = 18,  alpha = 0.6,
                   col =  beak_col)  + ggtitle(titl) + 
        theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
              plot.title = element_text(size = 17))+
        theme(legend.position = "none") + 
        scale_x_discrete(labels=c("All" = expression(All%->%Cur),
                                  "IUCN" = expression(IUCN%->%Cur),
                                  "Fut" = expression(Cur%->%Fut)))
      
      #y-axis lims
      IT_YL <- range(select(obs, Alpha))
      IT_YL_U <- IT_YL[2] + (abs(diff(IT_YL)) * 0.1)
      IT_YL_L <- IT_YL[1] - (abs(diff(IT_YL)) * 0.1)                                                          
      
      jorisB <- ggplot(data = obs, aes(Type, Alpha)) + 
        theme_bw() + xlab("") + ylab("")  +
        geom_point(size = 3,
                   shape = 18, col = "#716aa3") + 
        geom_line(data=subset(obs, 
                              Type %in% c("All", "IUCN", "Cur")),
                  group = 1, col = "#716aa3") + 
        geom_line(data=subset(obs, 
                              Type %in% c("Cur","Fut")),
                  group = 1, linetype = 2,
                  col = "#716aa3") +
        ylim(IT_YL_L, IT_YL_U)
      
      joris <- list(list(jorisA, jorisB), obs, beak_ZP_mat)
    }#eo if hist
    
  }#eo if beak
  
  return(joris)
  
}

#########HYPERVOLUME ANALYSES######################

###Internal functions#####

#internal function to randomly assign n extinct species and remove
int_null <- function(dat, n){
  
  #random sample of n sp to go extinct
  n2 <- sample(1:nrow(dat), n, replace = F)
  #remove extinct sp
  dat2 <- dat[-n2,]
  return(dat2)
}

#internal function to fit the hypervolumes, convex hulls and trees
int_null3 <- function(methodz, mat_beak_dat, trait_null_dat = NULL, 
                      axes_F  = NULL, 
                      method.hv = NULL, svmg = NULL, kde.band = NULL,
                      SPP = NULL,
                      treez = NULL, use_picante = TRUE,
                      use_phylo = use_phylo){
  
  if (methodz == "hyper"){
    
    if (!identical(colnames(mat_beak_dat), rownames(trait_null_dat))){
      stop("Order error 3")
    }
    
    if (method.hv == "svm"){
      BAT_beak_null <- BAT::kernel.build(mat_beak_dat, 
                                         trait_null_dat[,axes_F], 
                                         method = method.hv, 
                                         svm.gamma = svmg,
                                         samples.per.point = SPP,
                                         abund = FALSE, axes = 0)
    } else if (method.hv == "gaussian"){
      BAT_beak_null <- BAT::kernel.build(mat_beak_dat, 
                                         trait_null_dat[,axes_F], 
                                         method = method.hv, 
                                         kde.bandwidth = kde.band,
                                         samples.per.point = SPP,
                                         abund = FALSE, axes = 0)
    } else{
      stop ("method.hv should be one of svm or gaussian")
    }
    
    return(BAT_beak_null)
    
  } else if (methodz == "convex"){
    if (!identical(colnames(mat_beak_dat), rownames(trait_null_dat))){
      stop("Order error 1")
    }
    byrd4 <-  BAT::hull.build(mat_beak_dat, 
                              trait_null_dat[,axes_F], axes = 0)
    BAT_beak_null <- BAT::hull.alpha(byrd4)
    return(BAT_beak_null)
    
  } else if (methodz == "tree"){
    if (!use_picante){
      BAT_beak_null <- BAT::alpha(mat_beak_dat, 
                                  tree = treez) %>% as.vector()
    } else {
      #using picante::pd for now as BAT takes much longer and uses more memory
      #if phylogeny used, include the root, but nj dendrogram = unrooted
      if(!is.logical(use_phylo)) stop("use_phylo must be logical")
      if (!use_phylo){
        BAT_beak_null <- picante::pd(mat_beak_dat, tree = treez,
                                     include.root = FALSE)[1,1]
      } else {
        BAT_beak_null <- picante::pd(mat_beak_dat, tree = treez,
                                     include.root = TRUE)[1,1]
      }
    }
    BAT_beak_null_dispersion <-  BAT::dispersion(mat_beak_dat, tree = treez,
                                                 abund = FALSE, 
                                                 relative = FALSE) %>% as.vector()
    return(list(BAT_beak_null, BAT_beak_null_dispersion))
  }
  
}


#########################################################################
######For each period, create a PA matrix and then run one of the
#####methods. If a null run, first randomly set N species extinct###########
######################################################################

#n4 = number of species to go extinct
#rn4 = rowname for PA matrix (e.g. H2C)
#repn = number of replicate runs for observed hyper run
#observed = TRUE (for observed values) or FALSE (fur null values)
int_null4 <- function(dat_ind4, n4 = NULL, rn4, repn = NULL,
                      method, axes_F = NULL, method.hv = NULL, svmg = NULL, 
                      BAND = NULL,
                      SPP = NULL, dendro = NULL, use_picante = TRUE, 
                      use_phylo = use_phylo,
                      observed){
  
  if (observed){
    if (is.numeric(n4)) stop("No extinctions for observed")
  } else {
    if (!is.numeric(n4)) stop("Need extinction number")
  }
  
  if (observed & method == "hyper"){
    if (!is.numeric(repn)) stop("repn must be numeric")
  }
  
  #If not observed then randomise the species extinctions
  if (!observed){
    dat_ind4_null <- int_null(dat_ind4, n =  n4)
  } else {
    dat_ind4_null <- dat_ind4
  }
  
  sp_beak_4_null <- dat_ind4_null$species
  mat_beak_4_null <- matrix(1, nrow = 1, ncol = length(sp_beak_4_null))
  colnames(mat_beak_4_null) <- sp_beak_4_null
  rownames(mat_beak_4_null) <- rn4
  
  if (method == "hyper"){
    
    if (observed){
      
      BAT_beak_null_4 <- replicate(repn, 
                                   {int_null3(methodz = method, 
                                              mat_beak_dat = mat_beak_4_null,
                                              trait_null_dat = dat_ind4_null, 
                                              axes_F = axes_F,
                                              method.hv = method.hv, 
                                              svmg = svmg, kde.band = BAND,
                                              SPP = SPP)}, simplify = FALSE)
      
    } else {
      BAT_beak_null_4 <- int_null3(methodz = method, mat_beak_dat = mat_beak_4_null,
                                   trait_null_dat = dat_ind4_null, axes_F = axes_F,
                                   method.hv = method.hv, svmg = svmg, kde.band = BAND,
                                   SPP = SPP)
    } #eo if observed
    
    
  } else if (method == "convex"){
    BAT_beak_null_4 <- int_null3(methodz = method, mat_beak_dat = mat_beak_4_null,
                                 trait_null_dat = dat_ind4_null, axes_F = axes_F)
    
  } else if (method == "tree"){
    
    BAT_beak_tree_4 <- int_null3(methodz = method, mat_beak_dat = mat_beak_4_null,
                                 treez = dendro, 
                                 use_picante = use_picante,
                                 use_phylo = use_phylo)
    BAT_beak_null_4 <- BAT_beak_tree_4[[1]]
    BAT_beak_null_4_dispersion <-  BAT_beak_tree_4[[2]]
    
  }
  
  if (method %in% c("hyper", "convex")){
    return(BAT_beak_null_4)
  } else if (method == "tree"){
    return(list(BAT_beak_null_4, BAT_beak_null_4_dispersion))
  }
}#eo function


#########################################################
#####
######################################################

int_null2 <- function(dat_ind, prehistoric_F, historic_F, current_F, 
                      fut_N,
                      axes_F, method, method.hv = NULL, svmg = NULL, 
                      BAND = NULL,
                      SPP = NULL, dendro = NULL,
                      analysis,
                      use_picante = TRUE,
                      use_phylo = use_phylo){
  ##extinction numbers
  Hist_N <- table(dat_ind$ExtinctionPeriod)["IUCN"] %>% as.vector()
  Pre_N <- table(dat_ind$ExtinctionPeriod)["All"] %>% as.vector()
  All_N <- Hist_N + Pre_N
  #number of future extinctions = fut_N
  
  #P2C
  BAT_beak_null_C <- int_null4(dat_ind4 = prehistoric_F, n4 = All_N, 
                               rn4 = "All", method = method, 
                               axes_F = axes_F,
                               method.hv = method.hv, 
                               svmg = svmg, BAND = BAND,
                               SPP = SPP, dendro = dendro, 
                               use_picante = use_picante,
                               use_phylo = use_phylo,
                               observed = FALSE)
  
  #H2C
  BAT_beak_null_HC <- int_null4(dat_ind4 = historic_F, n4 = Hist_N, 
                                rn4 = "IUCN", method = method, 
                                axes_F = axes_F,
                                method.hv = method.hv, 
                                svmg = svmg, BAND = BAND,
                                SPP = SPP, dendro = dendro, 
                                use_picante = use_picante,
                                use_phylo = use_phylo,
                                observed = FALSE)
  
  #C2F
  BAT_beak_null_CF <- int_null4(dat_ind4 = current_F, n4 = fut_N, 
                                rn4 = "Fut", method = method, 
                                axes_F = axes_F,
                                method.hv = method.hv, svmg = svmg, 
                                BAND = BAND,
                                SPP = SPP, dendro = dendro,
                                use_picante = use_picante,
                                use_phylo = use_phylo,
                                observed = FALSE)
  
  if (method == "tree"){
    
    BAT_beak_null_C_dispersion <- BAT_beak_null_C[[2]]
    BAT_beak_null_C <- BAT_beak_null_C[[1]]
    
    BAT_beak_null_HC_dispersion <- BAT_beak_null_HC[[2]]
    BAT_beak_null_HC <- BAT_beak_null_HC[[1]]
    
    BAT_beak_null_CF_dispersion <- BAT_beak_null_CF[[2]]
    BAT_beak_null_CF <- BAT_beak_null_CF[[1]]
  }#eo if tree
  
  if (method == "hyper"){
    if (analysis == "main"){
      beak_r <- matrix(NA, nrow = 3, ncol = 2)
      colnames(beak_r) <- c("Alpha", "Dispersion")
      
      beak_r[,2] <- c(as.vector(BAT::kernel.dispersion(BAT_beak_null_C, 
                                                       func = "divergence", frac = 1)),
                      as.vector(BAT::kernel.dispersion(BAT_beak_null_HC, 
                                                       func = "divergence", frac = 1)),
                      as.vector(BAT::kernel.dispersion(BAT_beak_null_CF, 
                                                       func = "divergence", frac = 1)))
      
      
    } else if (analysis == "beak"){
      beak_r <- matrix(NA, nrow = 3, ncol = 1)
      colnames(beak_r) <- c("Alpha")
      
    }
    
    #Alpha values in both main and beak
    beak_r[,1] <- c(as.vector(BAT::kernel.alpha(BAT_beak_null_C)),
                    as.vector(BAT::kernel.alpha(BAT_beak_null_HC)),
                    as.vector(BAT::kernel.alpha(BAT_beak_null_CF)))
    
    
  } else if (method == "convex"){
    beak_r <- matrix(NA, nrow = 3, ncol = 2)
    colnames(beak_r) <- c("Alpha", "Dispersion")
    beak_r[,1] <- c(BAT_beak_null_C, 
                    BAT_beak_null_HC, BAT_beak_null_CF)
  } else if (method == "tree"){
    beak_r <- matrix(NA, nrow = 3, ncol = 2)
    colnames(beak_r) <- c("Alpha", "Dispersion")
    beak_r[,1] <- c(BAT_beak_null_C, 
                    BAT_beak_null_HC, BAT_beak_null_CF)
    beak_r[,2] <- c(BAT_beak_null_C_dispersion, 
                    BAT_beak_null_HC_dispersion, BAT_beak_null_CF_dispersion)
  }
  
  rownames(beak_r) <-  c("All", "IUCN", "Fut")
  beak_r <- as.data.frame(beak_r)
  return(round(beak_r, 3))
  
}#eo int_null2


######calculate Z and p-values for each based on output of null hyper
#plot_null = whether or not to plot the histograms of the null values
#for each test
null_hyper_ZP <- function(beak_r_obs, hurricane, plot_null = FALSE){
  
  resZP_beak <- matrix(nrow = 3*ncol(beak_r_obs), ncol = 6)
  colnames(resZP_beak) <- c("Z", "P", "Metric", "Period", "SES", "SES_P")
  k <- 1
  #set up new version of observed values with C as first row 
  #and then another C and then F
  ZP_obs <- rbind(beak_r_obs[c("Cur"), c("Alpha", "Dispersion")], 
                  beak_r_obs[c("Cur", "Fut"), c("Alpha", "Dispersion")])
  
  rownames(ZP_obs) <- c("All", "IUCN", "Fut")
  
  #if plot_null, create a list to save the ggplots
  if (plot_null){
    plot_null_list <- vector("list", 
                             length = ncol(beak_r_obs) * 3)
  }
  #double for loop to iterate across the two metrics and three time 
  #periods, extract
  #the relevant null values and observed value, and then run zP
  for (i in 1:3){
    M2 <- c("All", "IUCN", "Fut")[i]
    for (j in 1:ncol(beak_r_obs)){
      M <- c("Alpha", "Dispersion")[j]
      dum_ZP <- vapply(hurricane, function(x) x[M2, M], 
                       FUN.VALUE = numeric(1))
      #plot histogram of the null values
      if (plot_null){
        ggtitlz <- paste0(M2, "_", M)
        plot_null_list[[k]] <- ggplot() +
          geom_histogram(data = as.data.frame(dum_ZP), aes(dum_ZP),
                         bins = 25) + ggtitle(ggtitlz)
      }
      ob_ZP <- ZP_obs[M2, M]
      resZP_beak[k, 1:2] <- round(zP(dum_ZP, ob_ZP, type = "ES"), 3)
      resZP_beak[k, 3] <- M
      resZP_beak[k, 4] <- M2
      resZP_beak[k, 5:6] <- round(zP(dum_ZP, ob_ZP, type = "SES"), 3)
      k <- k +1
    }#eo j
  }#eo i
  
  if (plot_null){
    jpeg(file = "null_histograms.jpeg", width = 30, 
         height = (length(plot_null_list) / 3) * 10, 
         units = "cm", res = 300)
    gridExtra::grid.arrange(grobs = plot_null_list, 
                            nrow = ceiling(length(plot_null_list) / 3))
    dev.off()
  }
  resZP_beak <- as.data.frame(resZP_beak)
  return(resZP_beak)
}


####MAIN FUNCTION######### 

#allSp2_F = main allSp2 object
##future_F = one run of future extinctions from the
#main future dataset (from iucn_sim)
#GEOG = "All", "Yes" (for island endemics)
#analysis = either "main" or "beak"
#method = "hyper" or "convex" or "tree"
#observed = TRUE (observed metrics) or FALSE (null model runs)
#dendro = NULL or the global species dendrogram
# rn = number of replications to do for the observed hypervolume values (average then taken)

null_hyper <- function(allSp2_F = allSp2, 
                       future_F,
                       GEOG, analysis = "main",
                       axes_N = NULL,
                       method = "hyper", observed = TRUE,
                       method.hv = NULL, svmg = NULL, 
                       BAND = NULL,
                       SPP = NULL,
                       dendro = NULL,
                       use_picante = TRUE, use_phylo = FALSE,
                       rn = NULL){
  
  #select the correct PCA axes depending on the analysis
  if (analysis == "main"){
    axes_F <- c("PC1", "PC2", "PC3", "PC4", "PC5")
    axes_F <- axes_F[1:axes_N]
    #  axes_F <- c("PC1", "PC2", "PC3")
  } else if (analysis == "beak"){
    axes_F <- c("PC1_beak", "PC2_beak", "PC3_beak")
  }
  
  #this then filters or sets the dataset to the correct version based on selected
  #GEOG value
  if (GEOG == "All"){
    dat_ind <- allSp2_F
    prehistoric_F <- allSp2_F
    historic_F <- allSp2_F %>% 
      filter(ExtinctionPeriod %in% c("IUCN", "Current")) 
    current_F <- allSp2_F %>% 
      filter(ExtinctionPeriod == "Current")
    fut_N <- length(future_F[[1]])#all elements have same length
  } else {
    dat_ind <- filter(allSp2_F, IslandEndemic == GEOG)
    prehistoric_F <- allSp2_F %>% 
      filter(IslandEndemic == GEOG) 
    historic_F  <- allSp2_F %>% 
      filter(ExtinctionPeriod %in% c("IUCN", "Current") &
               IslandEndemic == GEOG)
    current_F  <- allSp2_F %>% 
      filter(ExtinctionPeriod == "Current" &
               IslandEndemic == GEOG)
    fut_N <- length(future_F[[1]])#all elements have same length
  }#eo if GEOG == ALL
  
  #########################
  ##Observed values
  ################################
  
  ##PA matrix for prehistoric period
  rownames(prehistoric_F) <- prehistoric_F$species
  rownames(historic_F) <- historic_F$species
  rownames(current_F) <- current_F$species
  
  if (observed){
    
    ##Prehistoric
    BAT_beak_obs_P <- int_null4(dat_ind4 = prehistoric_F, rn4 = "All", 
                                repn = rn, method = method, 
                                axes_F = axes_F, method.hv = method.hv, 
                                svmg = svmg, BAND = BAND,
                                SPP = SPP,
                                dendro = dendro, 
                                use_picante = use_picante, use_phylo = use_phylo,
                                observed = TRUE)
    
    ##historic
    BAT_beak_obs_H <- int_null4(dat_ind4 = historic_F, rn4 = "IUCN", 
                                repn = rn, method = method, 
                                axes_F = axes_F, 
                                method.hv = method.hv, 
                                svmg = svmg, BAND = BAND,
                                SPP = SPP, dendro = dendro,
                                use_picante = use_picante,
                                use_phylo = use_phylo,
                                observed = TRUE)
    
    ##Current
    BAT_beak_obs_C <- int_null4(dat_ind4 = current_F, rn4 = "Cur", 
                                repn = rn, method = method, 
                                axes_F = axes_F, 
                                method.hv = method.hv, 
                                svmg = svmg, BAND = BAND,
                                SPP = SPP, dendro = dendro,
                                use_picante = use_picante,
                                use_phylo = use_phylo,
                                observed = TRUE)
    
    ##Future
    
    #create new GEOG object to allow easy filtering
    if (GEOG == "All"){
      GEOG2 <- c("Yes", "No")
    } else if (GEOG == "Yes"){
      GEOG2 <- "Yes"
    } else {
      stop("Dualz")
    }
    
    if (!all(unlist(future_F) %in% 
             filter(allSp2_F, ExtinctionPeriod == "Current" &
                    IslandEndemic %in% GEOG2)$species)){
      stop("cigs")
    }
    #a list where each element is a list corresponding to
    #an individual run from future_F (iucn_sim runs). Each
    #of these lists holds the rn replicate hypervolumes etc.
    BAT_beak_obs_F <-  lapply(future_F, function(x){
      future1_F <- allSp2_F %>% 
        filter(ExtinctionPeriod == "Current" &
                 IslandEndemic %in% GEOG2 & (!species %in% x))
      if (!identical(rownames(future1_F),future1_F$species)){
        stop("grazze")
      }
      #rn set to 1 for these future runs due to computational demands.
      int_null4(dat_ind4 = future1_F, rn4 = "Fut", 
                repn = 1, method = method, 
                axes_F = axes_F, 
                method.hv = method.hv, 
                svmg = svmg, BAND = BAND,
                SPP = SPP, dendro = dendro,
                use_picante = use_picante,
                use_phylo = use_phylo,
                observed = TRUE)
    })
    
    
    ################################### 
    
    if (method == "hyper"){
      if (analysis == "main"){
        
        beak_r_obs <- matrix(NA, nrow = 4, ncol = 2)
        colnames(beak_r_obs) <- c("Alpha", "Dispersion")
        
        #take the mean value from the 5 replication runs
        hyper_disp_P <- vapply(BAT_beak_obs_P, 
                               function(x) as.vector(BAT::kernel.dispersion(x, func = "divergence", frac = 1)), 
                               FUN.VALUE = numeric(1)) %>% median()
        hyper_disp_H <- vapply(BAT_beak_obs_H, 
                               function(x) as.vector(BAT::kernel.dispersion(x, func = "divergence", frac = 1)), 
                               FUN.VALUE = numeric(1)) %>% median()
        hyper_disp_C <- vapply(BAT_beak_obs_C, 
                               function(x) as.vector(BAT::kernel.dispersion(x, func = "divergence", frac = 1)), 
                               FUN.VALUE = numeric(1)) %>% median()
        
        #this returns the median (across rn replicates) for each of
        #the iucn_sim runs; if only using 1 rn replicate, just returns the value
        hyper_disp_F1 <- lapply(BAT_beak_obs_F, function(y){
          vapply(y, function(x) as.vector(BAT::kernel.dispersion(x, func = "divergence", frac = 1)), 
                 FUN.VALUE = numeric(1)) %>% median()
        })
        #We then take the median of this value
        hyper_disp_F <- median(unlist(hyper_disp_F1))
        
        beak_r_obs[,2] <- c(hyper_disp_P, hyper_disp_H, 
                            hyper_disp_C, hyper_disp_F)
        
      } else if (analysis == "beak"){
        beak_r_obs <- matrix(NA, nrow = 4, ncol = 1)
        colnames(beak_r_obs) <- c("Alpha")
        
      }
      #Alpha hyper values (same code for main and beak analyses); take median alpha of N replications
      hyper_alpha_P <- vapply(BAT_beak_obs_P, 
                              function(x) as.vector(BAT::kernel.alpha(x)),
                              FUN.VALUE = numeric(1)) %>% median()
      hyper_alpha_H <- vapply(BAT_beak_obs_H, 
                              function(x) as.vector(BAT::kernel.alpha(x)), 
                              FUN.VALUE = numeric(1)) %>% median()
      hyper_alpha_C <- vapply(BAT_beak_obs_C, 
                              function(x) as.vector(BAT::kernel.alpha(x)), 
                              FUN.VALUE = numeric(1)) %>% median()
      #this returns the median (across rn replicates) for each of
      #the iucn_sim runs
      hyper_alpha_F1 <- lapply(BAT_beak_obs_F, function(y){
        vapply(y, function(x) as.vector(BAT::kernel.alpha(x)), 
               FUN.VALUE = numeric(1)) %>% median()
      })
      #We then take the median of this value
      hyper_alpha_F <- median(unlist(hyper_alpha_F1))
      
      beak_r_obs[,1] <- c(hyper_alpha_P, hyper_alpha_H, 
                          hyper_alpha_C, hyper_alpha_F)
      
      if (analysis == "main"){
        horsfall <- list("alpha" = hyper_alpha_F1,
                         "dispersion" = hyper_disp_F1)
      }
      
    } else if (method == "tree") {
      
      P_dispersion <- BAT_beak_obs_P[[2]]
      P_alpha <- BAT_beak_obs_P[[1]]
      
      H_dispersion <- BAT_beak_obs_H[[2]]
      H_alpha <- BAT_beak_obs_H[[1]]
      
      C_dispersion <- BAT_beak_obs_C[[2]]
      C_alpha <- BAT_beak_obs_C[[1]]
      
      #this returns the values for each of
      #the iucn_sim runs
      BAT_beak_obs_F1 <- vapply(BAT_beak_obs_F, function(x){
        x[[1]]
      }, FUN.VALUE = numeric(1))
      #We then take the median of this value
      F_alpha <- median(BAT_beak_obs_F1)
      #and same for dispersion
      F1_dispersion <- vapply(BAT_beak_obs_F, function(x){
        x[[2]]
      }, FUN.VALUE = numeric(1))
      F_dispersion <- median(F1_dispersion)
      
      beak_r_obs <- matrix(NA, nrow = 4, ncol = 2)
      colnames(beak_r_obs) <- c("Alpha", "Dispersion")
      beak_r_obs[,1] <- c(P_alpha, H_alpha, 
                          C_alpha, F_alpha)
      beak_r_obs[,2] <- c(P_dispersion, H_dispersion, 
                          C_dispersion, F_dispersion)
      
      horsfall <- list("alpha" = BAT_beak_obs_F1,
                       "dispersion" = F1_dispersion)
      
    } else if (method == "convex"){
      
      BAT_beak_obs_F_M <- median(unlist(BAT_beak_obs_F))
      beak_r_obs <- matrix(NA, nrow = 4, ncol = 2)
      colnames(beak_r_obs) <- c("Alpha", "Dispersion")
      beak_r_obs[,1] <- c(BAT_beak_obs_P, BAT_beak_obs_H, 
                          BAT_beak_obs_C, BAT_beak_obs_F_M)
      horsfall <- list("alpha" = BAT_beak_obs_F)
    }
    
    rownames(beak_r_obs) <-  c("All", "IUCN", "Cur", "Fut")
    beak_r_obs <- as.data.frame(beak_r_obs)
    
    if (analysis == "main") {
      return(list(beak_r_obs, horsfall))
    } else{
      return(beak_r_obs)
    }
    
  } else if (!observed){
    
    #########################
    ##Null values
    ################################
    
    hurricane <- int_null2(dat_ind = dat_ind, 
                           prehistoric_F = prehistoric_F, 
                           historic_F = historic_F, 
                           current_F = current_F,
                           fut_N = fut_N,
                           axes_F = axes_F, method = method, 
                           method.hv = method.hv,
                           svmg = svmg, BAND = BAND, 
                           SPP = SPP, dendro = dendro, 
                           analysis = analysis,
                           use_picante = use_picante, 
                           use_phylo = use_phylo)

    return(hurricane)
  }
  
}

#######################################################################
##function to plot null model results
#side = one or two sided P-test
# beak_run_F = main_res_list
# method = "hyper"
# side = "one" #one or two sided p-value
# BW = "nrd" #bandwidth for geom_violin 
# inset_only = FALSE #only plot the inset plot (only works for method = "hypervolume")
##########################################################

#colour pallete (colour-blind friendly):
#https://grafify-vignettes.netlify.app/colour_palettes.html

hyper_plot <- function(beak_run_F, ggtitl = "a) All ", 
                       method = "hyper",
                       side = "one", BW = "nrd", 
                       inset_only = FALSE){
  
  ##Extract and format the null values
  if (!inset_only){
    
    beak_list <- data.table::rbindlist(beak_run_F[[2]], 
                                       fill=FALSE, idcol=NULL)
    beak_list$Type <- rep(c("All", "IUCN", "Fut"), 
                          (nrow(beak_list) / 3))
    
    beak_list$Type <- factor(beak_list$Type, 
                             levels = c("All", "IUCN", "Fut"))
    
    if (method == "hyper" | method == "tree"){
      
      beak_long <- tidyr::pivot_longer(data = beak_list, 
                                       cols = c(Alpha, Dispersion),
                                       names_to = "Metric")
    } else if (method == "convex"){
      beak_long <- tidyr::pivot_longer(data = beak_list, 
                                       cols = c(Alpha), 
                                       names_to = "Metric")
    }
  }#eo if inset only
  
  #store original observed values for inset plot
  #if inset_only = TRUE, only the observed values are provided, and
  #so this is the dataframe to use, whereas otherwise the observed values
  #are the first element in the list
  if (inset_only){
    obs_three <- beak_run_F
  } else {
    obs_three <- beak_run_F[[1]]
  }
  
  obs_null2 <- mutate(obs_three, Type = c("All", "IUCN", "Cur", "Fut")) 
  obs_null2$Type <- factor(obs_null2$Type, 
                           levels = c("All", "IUCN", "Cur", "Fut"))
  if (method == "hyper" | method == "tree"){
    obs_null_2_long <- tidyr::pivot_longer(data = obs_null2, 
                                           cols = c(Alpha, Dispersion),
                                           names_to = "Metric")
  } else if (method == "convex"){
    obs_null_2_long <- tidyr::pivot_longer(data = obs_null2, 
                                           cols = c(Alpha), 
                                           names_to = "Metric")
  }
  
  ###y-axis limits
  #This is a long-winded approach as it needs to adjust the axis limits
  #in combination with facet_wraps (i.e. different limits for the 
  #different facets). It works by creating a new version of the observed
  #value dataframe, then replacing some of the observed values with the
  #respective axis limit values, and then using:
  #     geom_blank(data=ddummy, aes(Type, value)) 
  #in the plot to create a plot with the extended axis limits
  #https://stackoverflow.com/questions/30280499/different-y-limits-on-ggplot-facet-grid-bar-graph
  #2b) change inset only above
  
  ##alpha
  IT_YL_A <- obs_null_2_long %>%
    filter(Metric == "Alpha") %>%
    select(value) %>%
    range()
  IT_YL_U_A <- IT_YL_A[2] + (abs(diff(IT_YL_A)) * 0.1)
  IT_YL_L_A <- IT_YL_A[1] - (abs(diff(IT_YL_A)) * 0.1)
  ddummy <- obs_null_2_long
  ddummy$value[which(ddummy$Metric == "Alpha")][1:2] <-
    c(IT_YL_L_A, IT_YL_U_A)
  
  ##dispersion
  if (method == "hyper" | method == "tree"){
    IT_YL_D <- obs_null_2_long %>%
      filter(Metric == "Dispersion") %>%
      select(value) %>%
      range()
    IT_YL_U_D <- IT_YL_D[2] + (abs(diff(IT_YL_D)) * 0.1)
    IT_YL_L_D <- IT_YL_D[1] - (abs(diff(IT_YL_D)) * 0.1)
    ddummy$value[which(ddummy$Metric == "Dispersion")][1:2] <-
      c(IT_YL_L_D, IT_YL_U_D)
  }
  
  ##make small version just of the four observed values for inset plot
  joris2B <- obs_null_2_long %>% 
    ggplot(aes(Type, value)) + 
    geom_blank(data=ddummy, aes(Type, value)) +
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) + 
    theme(strip.background =element_rect(fill="white"),
          strip.text = element_text(face = "bold")) +
    xlab("") + ylab("")  +
    geom_point(size = 3,
               shape = 18, col = "#716aa3") +
    facet_wrap(~Metric, scales = "free") + 
    ggtitle(ggtitl)+ 
    geom_line(data=subset(obs_null_2_long, 
                          Type %in% c("All", "IUCN", "Cur")),
              group = 1, col = "#716aa3") + 
    geom_line(data=subset(obs_null_2_long, 
                          Type %in% c("Cur","Fut")),
              group = 1, linetype = 2,
              col = "#716aa3")
  
  if (inset_only){   
    
    return(joris2B)
    
  } #eo if inset_only
  
  #format the observed value (the first is the P value, 
  #so extract and remove and repeat the C value twice). 
  beak_obs_P <- beak_run_F[[1]][1,]
  slide_away <- cbind(rbind(beak_run_F[[1]][c("Cur"),
                                            c("Alpha", "Dispersion")], 
                            beak_run_F[[1]][c("Cur", "Fut"),
                                            c("Alpha", "Dispersion")]))
  
  if (method == "hyper" | method == "tree"){
    beak_run_F[[1]] <- slide_away
  } else if(method == "convex"){
    beak_run_F[[1]] <- slide_away[,"Alpha", drop = FALSE]
  } 
  
  beak_run_F[[1]]$Type <- c("All", "IUCN", "Fut")
  beak_run_F[[1]]$Type <- factor(beak_run_F[[1]]$Type, 
                                 levels = c("All", "IUCN", "Fut"))
  
  if (method == "hyper" | method == "tree"){
    beak_obs_long <- tidyr::pivot_longer(data = beak_run_F[[1]], 
                                         cols = c(Alpha, Dispersion), 
                                         names_to = "Metric")
  } else if (method == "convex"){
    beak_obs_long <- tidyr::pivot_longer(data = beak_run_F[[1]], 
                                         cols = c(Alpha), 
                                         names_to = "Metric")
  } 
  
  ##set point colour based on significance of ES
  if (method == "convex"){
    beak_run_F[[3]] <- filter(beak_run_F[[3]], 
                              Metric == "Alpha")
  }
  if (side == "two"){
    beak_obs_long$SIG <- ifelse(beak_run_F[[3]]$P > 0.975 | 
                                  beak_run_F[[3]]$P < 0.025,
                                "Sig", "NS")
  } else if (side == "one"){
    beak_obs_long$SIG <- ifelse(beak_run_F[[3]]$P < 0.05, 
                                "Sig", "NS")
  }
  
  #check everything in the right order
  if(!((identical(beak_run_F[[3]]$Metric, 
                  beak_obs_long$Metric)) &
       (identical(beak_run_F[[3]]$Period,
                  as.character(beak_obs_long$Type))))){
    stop("RS")
  }
  
  #set bandwith based on value of BW
  if (BW != "nrd") BW <- as.numeric(BW)
  
  #set y-axis title
  if (method == "tree"){
    ylab_div <- "Phylogenetic diversity"
  } else {
    ylab_div <- "Functional diversity"
  }
  
  ##violin colour
  # if (grepl("phylo", ggtitl)){
  #   vio_col <- "#F3D54D"
  #   vio_col <- cls2[2]
  # } else{
  #   vio_col <- "#48ABAE"
  vio_col <- "#E2E2E2"
  # }
  
  joris2A <- ggplot(data = beak_long, aes(Type, value)) + 
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) + 
    theme(strip.background =element_rect(fill="white"),
          strip.text = element_text(face = "bold")) +
    theme(axis.text = element_text(size=10),
          axis.title = element_text(size = 13)) +
    xlab("Time Period") + 
    ylab(ylab_div) + 
    geom_violin(aes(fill = Type), bw = BW,
                colour = vio_col)  +
    geom_point(data = beak_obs_long, aes(Type, value, 
                                         col = SIG), 
               size = 3.25,
               shape = 18,
               alpha = 0.8)  + ggtitle(ggtitl) +
    facet_wrap(~Metric, scales = "free") +
    #    scale_fill_grafify(palette = "fishy") +
    scale_fill_manual(values = rep(vio_col,4),
                      name = "Metric") +
    scale_color_manual(values = c("Sig" = "#023FA5", 
                                  "NS" =  "#8E063B"),
                       name = "Metric") +
    theme(legend.position = "none") + 
    scale_x_discrete(labels=c("All" = expression(All%->%Cur),
                              "IUCN" = expression(IUCN%->%Cur),
                              "Fut" = expression(Cur%->%Fut)))
  
  joris2 <- list(joris2A, joris2B)
  
  
  return(joris2)
}

################################################################
###############CONTRIBUTION FUNCTIONS##########################
#############################################################

#fit the hypervolume and return the contribution values for all species
#hyp_method = gaussian or svm
# SPP = samples per point argument for SVM hypervolumes
# rel = TRUE or FALSE (return relative contribution values
#that sum to 1 or not)
#axes_N = number of PCA axes to use
# FD_PD = whether to use FD and hypervolumes (FD_hyper), or FD/PD and trees
# tree = if FD_PD = "FD" or "PD", the phylogeny (or dendrogram)
# type = "ES" or "SES"
# null_contr_N = number of null iterations
contribution_own <- function(allSp_data = allSp2, 
                             hyp_meth = "svm",
                             SPP,
                             rel = FALSE, axes_N = 5, 
                             FD_PD, tree = NULL, 
                             type = "ES",
                             null_contr_N = 999){
  
  #PA matrix with all species (extant and extinct)
  paCont = matrix(1, ncol = nrow(allSp_data), nrow = 1)
  colnames(paCont) = rownames(allSp_data)
  
  if (FD_PD == "FD_hyper"){
    axes_F <- c("PC1", "PC2", "PC3", "PC4", "PC5")
    axes_F <- axes_F[1:axes_N]
    
    
    #Build the hypervolume with n PCA axes; is using all
    #species so don't need to set bandwidth and samples.per.point
    if (!identical(colnames(paCont), rownames(allSp_data))){
      stop("Order error 2")
    }
    kCont <- BAT::kernel.build(paCont, 
                               allSp_data[,axes_F], 
                               abund = FALSE, axes = 0,
                               method.hv = hyp_meth,
                               samples.per.point = SPP)
    #  method.hv = "svm")
    #using relative returns the contributions as a % of total
    #volume.
    k2Cont <- BAT::kernel.contribution(kCont, relative = rel)
    
    #add contributions as a column to allSp2
    if (!identical(names(k2Cont), rownames(allSp_data))) stop("Anjuna")
    allSp_data$contributions <-  k2Cont
    
    cont_p4 <- contr_tab(allSp_data, FD_PD = FD_PD, rel = rel)
    
    ##null model results
    #type = "ES" or "SES"
    nmrc4 <- null_model_results_contr(allSp_data, FD_PD = FD_PD,
                                      cont_p4, type = type,
                                      null_contr_N = null_contr_N,
                                      rel = rel)
    
    k3Cont <- list(cont_p4,  nmrc4,  k2Cont)
    
  } else if (FD_PD == "FD"| FD_PD == "PD"){
    #currently there is a bug in BAT which throws an error unless there
    #are more than one rows in the PA matrix
    paCont2 <- rbind(paCont, rep(0, ncol(paCont)))
    
    #for island endemics, prune tree to just have island endemics
    #We have to prune, as the alternative is to use the full tree,
    #and simply subset the island endemic species values. But this
    #means that the summed contribution values wont = total FD.
    
    if (ncol(paCont2) < 3500){
      tree2 <- ape::keep.tip(tree, colnames(paCont2))
    } else {
      tree2 <- tree
    }
    k2Cont <- BAT::contribution(paCont2, tree2, 
                                relative = rel)
    
    k2Cont <- k2Cont[1,]
    
    #This needs doing as phylo allSp2 was 
    #made using phylo_all[[1]] and the other phylos have same
    #tips but some in different places, and BAT::contribution
    #returns the values in order of tip labels in tree, and not
    #in order of column names in pa matrix!
    k2Cont <- k2Cont[match(colnames(paCont2), names(k2Cont))]
    
    #add contributions as a column to allSp2
    if (!identical(names(k2Cont), 
                   rownames(allSp_data))) stop("Anjuna")
    allSp_data$contributions <-  k2Cont
    
    cont_p4 <- contr_tab(allSp_data, 
                         FD_PD = FD_PD, rel = rel)
    ##null model results
    #type = "ES" or "SES"
    nmrc4 <- null_model_results_contr(allSp2_contr = allSp_data, 
                                      FD_PD = FD_PD,
                                      cont_p4, type = type,
                                      null_contr_N = null_contr_N,
                                      rel = rel)
    
    k3Cont <- list(cont_p4,  nmrc4,  k2Cont)
  }
  return(k3Cont)
}


#internal function
#extract TH (CR + EN + VU) and non threatened (NT),
#add to cont_P1 and remove the current row
cont_ext <- function(cont_dat = cont_P1B, cont_main = cont_P1A){
  P1b_TH <- filter(cont_dat, status %in% c("CR", "EN", "VU"))
  P1b_TH <- tibble("ExtinctionType" = "TH", "n" = sum(P1b_TH$n))
  P1b_NT <- filter(cont_dat, status %in% c("DD", "LC", "NT"))
  P1b_NT <- tibble("ExtinctionType" = "NT", "n" = sum(P1b_NT$n))
  cont_main <- rbind.data.frame(cont_main, P1b_TH, P1b_NT)
  cont_main <- cont_main[-which(cont_main$ExtinctionType == "Current"),]
  return(cont_main)
}


##build the table with % contribution for each category
#(EP, EPU, EH, CR, TH, NT)
contr_tab <- function(allSp_data, FD_PD, rel = FALSE){
  
  #plot dataset
  cont_P1A <- allSp_data %>%
    group_by(ExtinctionType) %>%
    summarise(n = sum(contributions))
  
  ##split into CR and other threatened
  cont_P1B <- allSp_data %>%
    group_by(status) %>%
    summarise(n = sum(contributions))
  
  cont_P1 <- cont_ext(cont_P1B, cont_P1A)
  
  if (rel){
    if (round(sum(cont_P1$n),5) != 1) stop("Fairport")
  }
  
  ##species richness totals
  cont_P2A <- allSp_data %>%
    group_by(ExtinctionType) %>%
    tally()
  
  #split into CR and other threatened
  cont_P2B <- allSp_data %>%
    group_by(status) %>%
    tally()
  
  cont_P2 <- cont_ext(cont_P2B, cont_P2A)
  #%
  cont_P2$n / nrow(allSp_data)
  
  cont_p3 <- rbind.data.frame(cont_P1, cont_P2)
  cont_p3$Contribution <- c(rep(FD_PD, 5), rep("SR", 5))

  cont_p3$ExtinctionType <- factor(cont_p3$ExtinctionType,
                                   levels = c("NT", "TH", 
                                              "IUCN", "EPA", "EPU"))
  
  ##second version with extinct species grouped
  cont_p4_A <- filter(cont_p3, 
                      ExtinctionType %in% c("IUCN", "EPA", "EPU"))
  
  cont_p4_B <- filter(cont_p3, 
                      ExtinctionType %in% c("TH", "NT"))
  
  cont_p4_sum <- cont_p4_A %>%
    group_by(Contribution) %>%
    summarise("Sumz" = sum(n))
  
  cont_p4_C <- data.frame("ExtinctionType" = rep("EX", 2),
                          n = as.vector(cont_p4_sum$Sumz),
                          Contribution = c(FD_PD, "SR"))
  
  cont_p4_D <- bind_rows(slice(cont_p4_C, 
                               which(Contribution == FD_PD)),
                         slice(cont_p4_B, 
                               which(Contribution == FD_PD)),
                         slice(cont_p4_C, 
                               which(Contribution == "SR")),
                         slice(cont_p4_B, 
                               which(Contribution == "SR")))
  
  if (cont_p4_D$n[1] != sum(cont_p3$n[1:3]) |
      cont_p4_D$n[4] != sum(cont_p3$n[6:8])){
    stop("Disco")
  }
  
  return(list(cont_p3, cont_p4_D))
}


##Null model contribution function
contr_null <- function(allSp_data, FD_PD, rel){
  
  #make a copy of allSp2, randomise the contribution values
  #across species and then calculate the % contributions across
  #each group (EP, CR etc)
  dum_allSp <- allSp_data
  dum_rn <- sample(dum_allSp$contributions, nrow(dum_allSp), 
                   replace = FALSE)
  dum_allSp$contributions <- dum_rn
  dum_contr <- contr_tab(dum_allSp, FD_PD, rel = rel)
  dum_res2 <- dum_contr[[2]][1:3,]$n
  names(dum_res2) <- dum_contr[[2]][1:3,]$ExtinctionType
  return(dum_res2)
}


##run the contribution null model and return the ZP table.
#type = "ES" or "SES"
null_model_results_contr <- function(allSp2_contr, FD_PD,
                                     cont_p4, type = "ES",
                                     null_contr_N = 999,
                                     rel){
  
  contr_rep <- replicate(null_contr_N, {contr_null(allSp2_contr, 
                                                   FD_PD, rel = rel)}, 
                         simplify = "matrix")
  
  contr_obs <- cont_p4[[2]]$n[1:3]#just the FD/PD results
  names(contr_obs) <- cont_p4[[2]]$ExtinctionType[1:3]
  
  if (!identical(names(contr_obs), rownames(contr_rep))) stop("Anjuna2")
  
  contr_zp <- matrix(ncol = 2, nrow = length(contr_obs))
  colnames(contr_zp) <- c("z", "p")
  rownames(contr_zp) <- names(contr_obs)
  
  for (i in 1:nrow(contr_rep)){
    contr_zp[i,] <- zP(dis = contr_rep[i,], 
                       obs = contr_obs[i], type = "ES")
  }
  contr_zp %>% round(3)
}

##Take the contribution output (either the 5 FD runs or the
#ten PD runs) and return the median and range of contribution
#values across these 5 runs, for each of the groups
null_format <- function(k5Cont){
  
  cont_p4 <- lapply(k5Cont, function(x) x[[1]][[2]])
  null_p4 <- lapply(k5Cont, function(x) x[[2]])
  
  cont_p4_M <- vapply(cont_p4, 
                      function(x) x$n[1:3], 
                      FUN.VALUE = numeric(3))
  
  cp4Med <- apply(cont_p4_M , 1, median) %>% round(2)
  cp4Ran <- apply(cont_p4_M , 1, range) %>% round(2)
  
  #relative values
  cont_p4_rel <- vapply(cont_p4, 
                        function(x) x$n[1:3] / 
                          sum(x$n[1:3]), 
                        FUN.VALUE = numeric(3))
  cp4Med_rel <- apply(cont_p4_rel , 1, median) %>% round(2)
  cp4Ran_rel <- apply(cont_p4_rel , 1, range) %>% round(2)
  
  cont_p4_T <- data.frame("ExtinctionType" = 
                            cont_p4[[1]]$ExtinctionType[1:3],
                          "Median" = cp4Med,
                          "Range_lower" = cp4Ran[1,],
                          "Range_upper" = cp4Ran[2,],
                          "Median_rel" =  cp4Med_rel,
                          "Range_lower_rel" = cp4Ran_rel[1,],
                          "Range_upper_rel" = cp4Ran_rel[2,])
  
  null_p4_MZ <- vapply(null_p4, 
                       function(x) as.data.frame(x)$z, 
                       FUN.VALUE = numeric(3))
  
  null_p4_MP <- vapply(null_p4, 
                       function(x) as.data.frame(x)$p, 
                       FUN.VALUE = numeric(3))
  
  np4MedZ <- apply(null_p4_MZ, 1, median) %>% round(2)
  np4RanZ <- apply(null_p4_MZ, 1, range) %>% round(2)
  np4MedP <- apply(null_p4_MP, 1, median) %>% round(2)
  np4RanP <- apply(null_p4_MP, 1, range) %>% round(2)
  
  null_p4_T <- data.frame("ExtinctionType" = rownames(as.data.frame(null_p4[[1]])),
                          "Median_z" = np4MedZ,
                          "Range_lower_z" = np4RanZ[1,],
                          "Range_upper_z" = np4RanZ[2,],
                          "Median_p" = np4MedP,
                          "Range_lower_p" = np4RanP[1,],
                          "Range_upper_p" = np4RanP[2,])
  return(list(cont_p4_T, null_p4_T))
}

##function to take the 50 pd contribution
#outputs and then look at how many times,
#for EX, TH and NT separately, the Z score is
#below zero (1) or above zero (2), and then
# if the P is significantly larger (1),
#significantly smaller (2) or non-sig (0)
contr_pd_table <- function(null_p4){
  soft <- do.call(rbind.data.frame, null_p4)
  soft$Type <- rep(c("EX", "TH", "NT"),
                   nrow(soft) / 3)
  hohne <- lapply(c("EX", "TH", "NT"), function(x){
    soft_f <- filter(soft, Type == x)
    landing <- ifelse(soft_f$z < 0, 1, 2)
    table(landing)
  })
  names(hohne) <- c("EX", "TH", "NT")
  
  hohne2 <- lapply(c("EX", "TH", "NT"), function(x){
    soft_f2 <- filter(soft, Type == x)
    landing2 <- sapply(soft_f2$p, function(x){
      if (x > 0.975){
        z2 <- 1
      } else if (x < 0.025){
        z2 <- 2
      } else {
        z2 <- 0
      }
      z2
    })
    table(landing2)
  })
  names(hohne2) <- c("EX", "TH", "NT")
  
  hohne3 <- list("Z" = hohne, "P" = hohne2)
  
  return(hohne3)
}

######################################################################
##function to calculate proportion of unique PD represented
#by extinct species in an order.
#sp = all (extant + extinct) species in an order
#sp2 = extinct species in an order
#tree = full global tree
#include_root = include.root argument in picante::pd
################################################################
PD_subset <- function(sp, sp2, tree, include_root = TRUE){
  
  ss <- matrix(1, nrow = 3, ncol = length(tree$tip.label))
  colnames(ss) <- tree$tip.label
  #1st row = all species
  ss[2,which(colnames(ss) %in% sp)] <- 0 #exclude all order species
  ss[3,which(colnames(ss) %in% sp2)] <- 0 #exclude just extinct order species
  
  ss_pd <- picante::pd(ss, tree, include.root = include_root)
  
  #unique branch length of order = total branch length -
  #(branch length of tree excluding all order species)
  ord_pd <- ss_pd$PD[1] - ss_pd$PD[2]
  
  #unique branch length of extinct order species =
  #total branch length - (branch length of tree excluding
  #extinct order species)
  ord_ext_pd <- ss_pd$PD[1] - ss_pd$PD[3] 
  
  #proportion of unique order PD represented by
  #unique extinct order PD
  ord_pd_res <- ord_ext_pd / ord_pd
  ord_pd_res <- round(ord_pd_res,2)
  
  return(ord_pd_res)
}


#########################################################
##############PD Results Formatting###########################
#############################################################

###Takes the PD output and formats for the main results
#table. Returns a list with two elements: (i) the median and
#range values for the percentage change and Z and P-values,
#for both alpha and dispersion, and (ii) the table significance
#results across the 50 trees:
#1 = sig lower than expected
#2 = sig higher than expected
#0 = non-sig
#Note that it also returns the median P2F values, for use in the
#temporal change plot
PD_multi_format <- function(main_run_1){
  
  ##observed values: return median / range perc_change
  pmtr1 <- lapply(main_run_1[[1]], function(x){
    pmtr1a <- perc_change(x$Alpha)
    pmtr1b <- perc_change(x$Dispersion)
    pmtr1c <- perc_change(c(x["All","Alpha"], 
                            x["Fut","Alpha"], 99, 99))[1]
    list(pmtr1a, pmtr1b, pmtr1c)
  })
  
  pmtr1_SS_A <- lapply(pmtr1, function(x) x[[1]]) %>%
    do.call(rbind.data.frame, .) %>%
    apply(., 2, function(y){
      c(median(y),
        range(y))}, simplify = "FALSE") %>%
    unlist() %>% as.vector()
  
  names(pmtr1_SS_A) <- c("PH_Alpha_Perc_Med", "PH_Alpha_Perc_RL",
                         "PH_Alpha_Perc_RU",
                         "PC_Alpha_Perc_Med", "PC_Alpha_Perc_RL",
                         "PC_Alpha_Perc_RU",
                         "HC_Alpha_Perc_Med", "HC_Alpha_Perc_RL",
                         "HC_Alpha_Perc_RU",
                         "CF_Alpha_Perc_Med", "CF_Alpha_Perc_RL",
                         "CF_Alpha_Perc_RU")
  
  pmtr1_SS_D <- lapply(pmtr1, function(x) x[[2]]) %>%
    do.call(rbind.data.frame, .) %>%
    apply(., 2, function(y){
      c(median(y),
        range(y))}, simplify = "FALSE") %>%
    unlist() %>% as.vector()
  
  names(pmtr1_SS_D) <- c("PH_Disp_Perc_Med", "PH_Disp_Perc_RL",
                         "PH_Disp_Perc_RU",
                         "PC_Disp_Perc_Med", "PC_Disp_Perc_RL",
                         "PC_Disp_Perc_RU",
                         "HC_Disp_Perc_Med", "HC_Disp_Perc_RL",
                         "HC_Disp_Perc_RU",
                         "CF_Disp_Perc_Med", "CF_Disp_Perc_RL",
                         "CF_Disp_Perc_RU")
  
  pmtr1_SS_P2F <- lapply(pmtr1, function(x) x[[3]]) %>%
    unlist() %>%
    median()
  
  names(pmtr1_SS_P2F) <- "PF_Alpha_Perc_Med"
  
  ##null model values: return median / range ES vals and p-vals
  #and also number of sig +- and non-sig cases for each
  pmtr2_A <- matrix(ncol = 6,
                    nrow =  length(main_run_1[[3]]))
  pmtr2_D <- matrix(ncol = 6,
                    nrow =  length(main_run_1[[3]]))
  
  for (i in 1:nrow(pmtr2_A)){
    dumPMA <- filter(main_run_1[[3]][[i]],
                     Metric == "Alpha")
    pmtr2_A[i,1:3] <- as.numeric(dumPMA$Z)
    pmtr2_A[i,4:6] <- as.numeric(dumPMA$P)
    
    dumPMD <- filter(main_run_1[[3]][[i]],
                     Metric == "Dispersion")
    pmtr2_D[i,1:3] <- as.numeric(dumPMD$Z)
    pmtr2_D[i,4:6] <- as.numeric(dumPMD$P)
  }
  
  pmtr2_SS_A <- pmtr2_A %>%
    apply(., 2, function(y){
      c(median(y),
        range(y))}, simplify = "FALSE") %>%
    unlist() %>% as.vector()
  
  pmtr2_SS_D <- pmtr2_D %>%
    apply(., 2, function(y){
      c(median(y),
        range(y))}, simplify = "FALSE") %>%
    unlist() %>% as.vector()
  
  names(pmtr2_SS_A) <- c("PC_Alpha_Z_Med", "PC_Alpha_Z_RL",
                         "PC_Alpha_Z_RU",
                         "HC_Alpha_Z_Med", "HC_Alpha_Z_RL",
                         "HC_Alpha_Z_RU",
                         "CF_Alpha_Z_Med", "CF_Alpha_Z_RL",
                         "CF_Alpha_Z_RU",
                         "PC_Alpha_P_Med", "PC_Alpha_P_RL",
                         "PC_Alpha_P_RU",
                         "HC_Alpha_P_Med", "HC_Alpha_P_RL",
                         "HC_Alpha_P_RU",
                         "CF_Alpha_P_Med", "CF_Alpha_P_RL",
                         "CF_Alpha_P_RU")
  
  names(pmtr2_SS_D) <- c("PC_Disp_Z_Med", "PC_Disp_Z_RL",
                         "PC_Disp_Z_RU",
                         "HC_Disp_Z_Med", "HC_Disp_Z_RL",
                         "HC_Disp_Z_RU",
                         "CF_Disp_Z_Med", "CF_Disp_Z_RL",
                         "CF_Disp_Z_RU",
                         "PC_Disp_P_Med", "PC_Disp_P_RL",
                         "PC_Disp_P_RU",
                         "HC_Disp_P_Med", "HC_Disp_P_RL",
                         "HC_Disp_P_RU",
                         "CF_Disp_P_Med", "CF_Disp_P_RL",
                         "CF_Disp_P_RU")
  
  #table of number of significant cases;
  #1 = sig lower, 2 = sig higher, 3 = non-sig
  pmtr3_A <- apply(pmtr2_A[,4:6], 2, function(y){
    soin1 <-  sapply(y, function(x){
      if (any(x < 0) | any(x > 1)) stop("post")
      if (x <0.05){
        z <- 1
      } else if (x > 0.95){
        z <- 2
      } else {
        z <- 0
      }
      z})
    soin1 <- factor(soin1, levels = 0:2)
    table(soin1)
  }, simplify = FALSE)
  
  if (!all(sapply(pmtr3_A, sum) == 
           nrow(pmtr2_A))) stop("ole")
  
  names(pmtr3_A) <- c("P2C_P", "H2C_P",
                      "C2F_P")
  
  ##Individual future simulation analyses
  fut_pd_res <- vector("list", length = 3)
  fut_pd_res[[1]] <- vector("list", length(main_run_1[[1]]))
  fut_pd_res[[2]] <- vector("list", length(main_run_1[[1]]))
  fut_pd_res[[3]] <- vector("list", length(main_run_1[[1]])) 
  #run across alpha and dispersion
  for (k in 1:2){
    #run across each phylogeny
    for (i in 1:length(main_run_1[[1]])){
      fut_pd_in <- main_run_1[[1]][[i]] #obs values for ith phylogeny
      #note there, the Fut value is the median F value across
      #the 100 simulations, using the ith phylogeny
      
      fut_pd_pc <- vector(length = length(main_run_1[[4]][[i]][[k]]))
      levP <- vector(length = length(main_run_1[[4]][[i]][[k]]))
      #Run across the 100 Future values
      for (j in 1:length(main_run_1[[4]][[i]][[k]])){
        #replace the F value in the obs table with the
        #jth individual F value
        fut_pd_in[,k][4] <- main_run_1[[4]][[i]][[k]][j]
        #Calculate and store the % change using this F-value
        fut_pd_pc[j] <- perc_change(fut_pd_in[,k])[4]
        
        if (k == 1){
          lev <- unlist(lapply(main_run_1[[2]][[i]], 
                               function(x){x[,k][3]}))
          levP[j] <- zP(lev, fut_pd_in[,k][4])[2]
        }
        
      }#eo j
      fut_pd_res[[k]][[i]] <- fut_pd_pc
      if (k ==1) fut_pd_res[[3]][[i]] <- levP
    }#eo i
  }#eo k
  
  ##Results vector / list
  pmtr4 <- c(pmtr1_SS_P2F, pmtr1_SS_A, pmtr1_SS_D, 
             pmtr2_SS_A, pmtr2_SS_D) 
  
  pmtr5 <- list(pmtr4, pmtr3_A, fut_pd_res)
  
  return(pmtr5)
}#eo function

#############################################################
##FIGURE 1 MAP
############################################################

world <- ne_countries(scale = "medium", returnclass = "sf") # Get a world map
world <- dplyr::filter(world, !region_un == "Antarctica") # Remove Antarctica

build_map <- function(extinct = extinct){
  
  
  #' For each area, find the number of species that went extinct either EP or EH
  #' We need a df with the area then the number of each extinction type, plus the lat/lon
  df1 <- filter(extinct, ExtinctionPeriodMap == "EX_Other") # Split the data into EH and EP
  df2 <- filter(extinct, ExtinctionPeriodMap %in% c("EX_IUCN",
                                                    "EX_POST1500"))
  ### Total species on each archipelago
  
  #combine the three sets of island group columns
  #this is because some species are on > 1 archipelago
  e1 <- dplyr::select(extinct, Isl_group_map:LON)
  e2 <- dplyr::select(extinct, Isl_group_map2:LON2) %>%
    rename(Isl_group_map = Isl_group_map2,
           LAT = LAT2, LON = LON2)
  e3 <- dplyr::select(extinct, Isl_group_map3:LON3) %>%
    rename(Isl_group_map = Isl_group_map3,
           LAT = LAT3, LON = LON3)
  eall <- rbind.data.frame(e1, e2, e3)
  all_N <- eall %>%
    group_by(Isl_group_map) %>%
    summarise("N" = n())
  #remove the NA row
  eall_NA <- which(is.na(all_N$Isl_group_map))
  if (length(eall_NA) !=1) stop("sunday morning")
  all_N <- all_N[-eall_NA,] 
  #add lat and lon to all_N for each unique allipelago
  all_M <- match(all_N$Isl_group_map, extinct$Isl_group_map)
  all_N$LAT <- extinct$LAT[all_M]
  all_N$LON <- extinct$LON[all_M]
  all_N$IslandEndemic <- extinct$IslandEndemic[all_M]
  
  ### Just Pre-1500
  
  #combine the three sets of island group columns
  #this is because some species are on > 1 archipelago
  e1 <- dplyr::select(df1, Isl_group_map:LON)
  e2 <- dplyr::select(df1, Isl_group_map2:LON2) %>%
    rename(Isl_group_map = Isl_group_map2,
           LAT = LAT2, LON = LON2)
  e3 <- dplyr::select(df1, Isl_group_map3:LON3) %>%
    rename(Isl_group_map = Isl_group_map3,
           LAT = LAT3, LON = LON3)
  eall <- rbind.data.frame(e1, e2, e3)
  df1_arch <- eall %>%
    group_by(Isl_group_map) %>%
    summarise("Pre1500" = n())
  #remove the NA row
  eall_NA <- which(is.na(df1_arch$Isl_group_map))
  if (length(eall_NA) !=1) stop("sunday morning")
  df1_arch <- df1_arch[-eall_NA,]
  #add lat and lon to df1_arch for each unique archipelago
  arch_M <- match(df1_arch$Isl_group_map, df1$Isl_group_map)
  df1_arch$LAT <- df1$LAT[arch_M]
  df1_arch$LON <- df1$LON[arch_M]
  df1_arch$IslandEndemic <- df1$IslandEndemic[arch_M]
  
  ### Just Post-1500
  
  #combine the three sets of island group columns
  #this is because some species are on > 1 archipelago
  e1 <- dplyr::select(df2, Isl_group_map:LON)
  e2 <- dplyr::select(df2, Isl_group_map2:LON2) %>%
    rename(Isl_group_map = Isl_group_map2,
           LAT = LAT2, LON = LON2)
  e3 <- dplyr::select(df2, Isl_group_map3:LON3) %>%
    rename(Isl_group_map = Isl_group_map3,
           LAT = LAT3, LON = LON3)
  eall <- rbind.data.frame(e1, e2, e3)
  df2_arch <- eall %>%
    group_by(Isl_group_map) %>%
    summarise("Post1500" = n())
  #remove the NA row
  eall_NA <- which(is.na(df2_arch$Isl_group_map))
  if (length(eall_NA) !=1) stop("sunday morning")
  df2_arch <- df2_arch[-eall_NA,]
  #add lat and lon to df2_arch for each unique archipelago
  arch_M <- match(df2_arch$Isl_group_map, df2$Isl_group_map)
  df2_arch$LAT <- df2$LAT[arch_M]
  df2_arch$LON <- df2$LON[arch_M]
  df2_arch$IslandEndemic <- df2$IslandEndemic[arch_M]
  
  arch_N <- merge(df1_arch, df2_arch, by = "Isl_group_map", all = T)# Join the datasets together
  arch_N <- arch_N[,c(1,2,6)] # Subset the island name, EP and EH cols
  arch_N[is.na(arch_N$Pre1500),]$Pre1500 <- 0  # Replace NA with 0
  arch_N[is.na(arch_N$Post1500),]$Post1500 <- 0 # Replace NA with 0
  arch_N <- dplyr::left_join(arch_N, all_N, 
                             by = "Isl_group_map") # Join the datasets
  arch_N$Pre1500_con <- 0                                                             
  arch_N$Post1500_con <- 0                                                             
  arch_N$Pre1500_con[15:18] <- arch_N$Pre1500[15:18]                                       
  arch_N$Post1500_con[15:18] <- arch_N$Post1500[15:18]                                       
  arch_N$Pre1500[15:18] <- 0                                                          
  arch_N$Post1500[15:18] <- 0                                                         
  sea_col <- rgb(140,218,255, maxColorValue = 255)# Sea colour             
  yel_col <- rgb(254,242,58, maxColorValue = 255)                               
  brow_col <- rgb(216,100,77, maxColorValue = 255)                               
  col1 <- "#3F0E14"
  col2 <- "#F2FFCC"
  
  arch_N[nrow(arch_N) + 1,] <- c("Legend1", 0, 1, 1, -55.299488, -8.047559, "No", 0,0)
  arch_N[nrow(arch_N) + 1,] <- c("Legend2", 0, 20, 20, -55.199291, 25.702440, "No", 0,0)
  arch_N[nrow(arch_N) + 1,] <- c("Legend3", 0, 80, 80, -59.373827, 2.850878, "No", 0,0)
  arch_N[,c(2:6,8,9)] <- apply(arch_N[,c(2:6,8,9)], 2, as.numeric)
  arch_N$N <- arch_N$N + 1 # Add one to help highlight small values
  
  
  
  g8 <- ggplot(data = world) + 
    theme_void() +
    geom_sf(color = "white", fill = "white") +
    theme(panel.background = element_rect(fill = sea_col),
          legend.position = "bottom",
          plot.title = element_blank()) +
    xlab("") + 
    ylab("") +
    geom_scatterpie(data = arch_N,aes(x=LON, y=LAT, r=2*log(N)),
                    cols = c("Pre1500", "Post1500", 
                             "Pre1500_con", "Post1500_con"), 
                    alpha = 0.7) +
    scale_fill_manual(values = c(brow_col, yel_col, col1, col2), 
                      guide = "none") +
    geom_scatterpie_legend(seq(1, ceiling(2*log(max(arch_N$N)))), 
                           x=-160, y=-55) 
  
  return(g8)
}
