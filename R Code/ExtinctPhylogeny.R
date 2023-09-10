
# Packages
library(treeio)
library(ggtree)
library(phytools)
require(parallel)
require(foreach)
library(snow)
library(doParallel)
library(doSNOW)
library(readxl)
library(MonoPhy)
library(ggplot2)
library(dplyr)
library(stringr)

##Function to graft a species to a given (truncated)
##branch
bind <- function(tree, node, per, per_fixed = FALSE, sp_name){
  #tree = tree object
  #node = target node
  #per = The fraction (0-1) of total branch length to truncate at
  #either end of the branch for grafting (e.g. 0.2 cuts
  #of 20% of the total branch lenth from either end and then
  #randomly selects the middle 60% for grafting)
  #per_fixed = whether to graft on at an exact place rather than random;
  #value between 0-1, with larger number meaning grafting happens closer
  #to the root.
  #sp_name = name of grafted species
  
  # Get the branch length
  Lx <- tree$edge.length[which(tree$edge[,2]==node)]   
  
  if (!per_fixed){
  #truncate the branch length
    LxTrun <- c((Lx * per), (Lx * (1 - per)))
  } else {
    LxTrun <- rep((Lx * per), 2)
  }

  # Bind the extinct sp.
  tree <- bind.tip(tree,                                                                   
                   paste0(sp_name), 
                   where = node, 
                   position = runif(1, min = LxTrun[1],
                                    max = LxTrun[2]))
  return(tree)
}


#Read in trees (This would be a set of BirdTree trees from the posterior distribution of BirdTree)
ctrees <- readRDS("Data/Jetz100.rds")
#subset a number you want to test
ctrees <- sample(ctrees, 50, replace = F)
length(ctrees[[1]]$tip.label) # 9993 tips (full tree)

## Add species in tree:
# We do it for N trees.
Ntree <- 50

## Read in Jetz master taxonomy (provided by BirdTree)
jetz <- read.csv("Data/Jetzmaster_taxonomy.csv")

##Set the percentage / fraction for branch truncation
#for the random grafting (see the above function for
#details)
PER <- 0.2

# Set number of cores
n.cores <- 20
cluster.ips <- NULL

if (is.null(cluster.ips)) {
  if (n.cores == 1) {
    `%dopar%` <- foreach::`%do%`
    on.exit(`%dopar%` <- foreach::`%dopar%`)
    cluster.ips <- NULL
  }
  else {
    temp.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  }
}else {
  Sys.setenv(R_PARALLEL_PORT = cluster.port)
  cluster.spec <- cluster_specification(cluster.ips = cluster.ips, 
                                        cluster.cores = cluster.cores, cluster.user = cluster.user)
  if (verbose == TRUE) {
    outfile <- ""
  }else {
    if (.Platform$OS.type == "windows") {
      outfile <- "nul:"
    }else {
      outfile <- "/dev/null"
    }
  }
  temp.cluster <- snow::makeSOCKcluster(ncores)
}
if (exists("temp.cluster")) {
  doParallel::registerDoParallel(cl = temp.cluster)
  doSNOW::registerDoSNOW(temp.cluster)
}

#Set a progress bar to return progress of the foreach loop
pb <- txtProgressBar(min = 0, max = length(Ntree), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

## Load in extinct species data
ex_data <- as.data.frame(readxl::read_xlsx("Data/Master_phylo_spreadsheet_V12.xlsx", sheet = 1))
ex_data$Group <- as.numeric(ex_data$Group)#NA warning fine (just because already NAs in Group)


#Run the parallel dataprep
ctreesComplete <- foreach(
  i = 1:Ntree,
  .options.snow = opts, 
  .packages = c("phytools")) %dopar%
  {
    
    ctree <- ctrees[[i]]   # Each loop we do one tree
  

    ## Reorder the dataset 
    ex <- ex_data
    ex <- ex[order(ex$phylo_id),]
    row.names(ex) <- 1:nrow(ex)
    
    ## Subset the species to randomly shuffle
    shuff <- dplyr::filter(ex, phylo_id2 == "xS")
    ## Remove those species from the initial DB
    ex <- dplyr::filter(ex, !phylo_id2 == "xS")
    ## Set the order
    groups <- as.numeric(unique(shuff$Group))
    groups <- sort(groups)
    
    for(i in 1:length(groups)){
      
      shuff2 <- dplyr::filter(shuff, Group == groups[i])
      shuff2 <-  shuff2[sample(1:nrow(shuff2)), ]
      ex <- rbind(ex, shuff2)
      
    }
    
    row.names(ex) <- 1:nrow(ex)
    
    # For each extinct species find the optimum place to bind 
    for(j in 1:nrow(ex)){
      
      ##' AP = Already present - Nothing needs to be done  
      if(ex$Type[j] == "AP"){
          next
      } else if(ex$Type[j] == "S"){
        ##' Scenario 1.1: Add target species as a single sister of species X ## - S (SISTER SPECIES)
        
        # Get the tip location for the sister sp.
        nodeX <- which(ctree$tip.label == paste0(ex$Sister_genus[j], "_", ex$Sister_species[j] ) ) 

        ##For the clades with many closely related extinct species,
        #set per to 0.75 to try to avoid really short terminal branches
        #(although sometimes this is forced due to BirdTree typology)
        if (ex$Group[j] %in% c(2, 4, 21, 25.1, 25.3, 34, 36)){
        
        # Bind the extinct sp.
        ctree <- bind(tree = ctree, node = nodeX,
                      per = 0.75, per_fixed = TRUE,
                      sp_name = ex$species[j])

        } else {
          # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])
        }

      } else if(ex$Type[j] %in% c("SSG", "SGG", "SGG2", "SFG", "SOG")){
        ##' Scenario 1.2: Add species as a sister (outgroup) of a group of species ## - SSG (SISTER SPECIES GROUP) & SGG 
        ##'               (SISTER GENUS GROUP) & SFG (SISTER FAMILY GROUP) & SOG (SISTER ORDER GROUP)
        
        if(ex$Type[j] == "SSG"){
          
          # Separate the species in the "sister_species_group" column
          sp <- stringr::str_split(ex$Sister_species_group[j], pattern = ";")
          spv <- vector()
          
          for(x in 1:length(sp[[1]])){
            
            spv <- c(spv, paste0( sp[[1]][x]) )
            
          } # Make the species group as a vector
          
          # This selects the most recent common ancestor for the group of species
          nodeX <- getMRCA(ctree, spv) 
          
          # Bind the extinct sp.
           ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])
          
        }else if(ex$Type[j] == "SGG"){
          
          # If only one species is present within the genus, then make a sister to that species
          if( length(ctree$tip.label[grep(paste0(ex$Sister_genus[j], "_"), ctree$tip.label)] ) == 1){
            
            # Get the tip location for the sister sp.
            nodeX <- which( ctree$tip.label == ctree$tip.label[grep(paste0(ex$Sister_genus[j], "_"), ctree$tip.label)] ) 
           
            # Bind the extinct sp.
            ctree <- bind(tree = ctree, node = nodeX,
                          per = PER, sp_name = ex$species[j])
          }else{
          
          # Get most recent common ancestor of genus 
          nodeX <- getMRCA( ctree, ctree$tip.label[grep(paste0(ex$Sister_genus[j], "_"), ctree$tip.label)] ) 
          
          # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])
          
          }
          
        }else if(ex$Type[j] == "SGG2"){
          
          sp <- stringr::str_split(ex$Sister_genus[j], pattern = ";")
          spv <- vector()
          
          # Get all the species in the genera
          for(x in 1:length(sp[[1]])){
            
            spv <- c(spv, ctree$tip.label[grep(paste0(sp[[1]][x], "_"), ctree$tip.label)])
            
          }
          
          # Get most recent common ancestor of species group
          nodeY <- getMRCA(ctree, spv)  
         
          # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeY,
                        per = PER, sp_name = ex$species[j])
          
        }else if(ex$Type[j] == "SFG"){
          
          # Get all species within the family 
          fam <- dplyr::filter(jetz, BLFamilyLatin == ex$Sister_family[j])
          spv <- vector()
          
          for(x in 1:nrow(fam)){
            
            spv <- c(spv, paste0(fam$Genus[x], "_", fam$Species[x]) )
            
          } # Make the species group as a vector
          
          # This selects the most recent common ancestor for the group of species
          nodeX <- getMRCA(ctree, spv)
         
          # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])

        }else if(ex$Type[j] == "SOG"){
          
          # Get all species within the family 
          ord <- dplyr::filter(jetz, Order == ex$Sister_order[j])
          spv <- vector()
          
          for(x in 1:nrow(ord)){
            
            spv <- c(spv, paste0(ord$Genus[x], "_", ord$Species[x]) )
            
          } # Make the species group as a vector
          
          # This selects the most recent common ancestor for the group of species
          nodeX <- getMRCA(ctree, spv) 
  
          # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])

        }
        
      } else if(ex$Type[j] %in% c("RSG", "RGG", "RGG2", "RCG", "RFG", "ROG")){
        ##' Scenario 2.1: Add species randomly within of a group of species ## - RSG (RANDOM SPECIES GROUP) & RGG 
        ##'               (RANDOM GENUS GROUP) & RFG (RANDOM FAMILY GROUP) & ROG (RANDOM ORDER GROUP)
        
        if(ex$Type[j] == "RSG"){
          
          # Separate the species in the "sister_species_group" column
          sp <- stringr::str_split(ex$Sister_species_group[j], pattern = ";")
          spv <- vector()
          
          for(x in 1:length(sp[[1]])){
            
            spv <- c(spv, paste0( sp[[1]][x]))
            
          } 
          
          ## Randomly select one of the species from the listed group. If it is not present in the tree yet, remove from the 
          ## list and select again. Will break if all species have been attempted and there were none present in the tree. 
          repeat{
            #Randomly select one of the species
            spv2 <- sample(spv, 1)
            # Get the tip location for the sister sp.
            nodeX <- which(ctree$tip.label == spv2) 
            # Break if there is a node value
            if(length(nodeX) != 0) break else{
              spv <- spv[!spv == spv2]
            }
            if(length(spv) == 0) break
            
          }
        
          ## Check if the node still is length zero
          if(length(nodeX) == 0){print(paste0("Node is still zero length for ", ex$Type[j], " for species ", ex$species[j],
                                                " (row ", j, ") after random species selection."))}
          
          ##For the clades with many closely related extinct species,
          #set per to 0.75 to try to avoid really short terminal branches
          #(although sometimes this is forced due to BirdTree typology)
          if (ex$Group[j] %in% c(5, 25.3, 53)){
            
            # Bind the extinct sp.
            ctree <- bind(tree = ctree, node = nodeX,
                          per = 0.75, per_fixed = TRUE,
                          sp_name = ex$species[j])
            
          } else {
            # Bind the extinct sp.
            ctree <- bind(tree = ctree, node = nodeX,
                          per = PER, sp_name = ex$species[j])
          }

        }else if(ex$Type[j] == "RGG"){
          
          # Get all the species in the genus
          sp <- stringr::str_split(ex$Sister_genus[j], pattern = ";")
          spv <- vector()
          
          # Get all the species in the genera
          for(x in 1:length(sp[[1]])){
            
            spv <- c(spv, ctree$tip.label[grep(paste0(sp[[1]][x], "_"), ctree$tip.label)])
            
          }
          
          #Randomly select one of the species
          spv2 <- sample(spv, 1)
          
          # Get the tip location for the sister sp.
          nodeX <- which( ctree$tip.label == spv2 ) 
          
          ##For the clades with many closely related extinct species,
          #set per to 0.75 to try to avoid really short terminal branches
          #(although sometimes this is forced due to BirdTree typology)
          if (ex$Group[j] %in% c(24, 25.2)){
            
            # Bind the extinct sp.
            ctree <- bind(tree = ctree, node = nodeX,
                          per = 0.75, per_fixed = TRUE,
                          sp_name = ex$species[j])
            
          } else {
            # Bind the extinct sp.
            ctree <- bind(tree = ctree, node = nodeX,
                          per = PER, sp_name = ex$species[j])
          }
          
        }else if(ex$Type[j] == "RGG2"){
          
          sp <- stringr::str_split(ex$Sister_genus[j], pattern = ";")
          spv <- vector()
          
          # Get all the species in the genera
          for(x in 1:length(sp[[1]])){
            
            spv <- c(spv, ctree$tip.label[grep(paste0(sp[[1]][x], "_"), ctree$tip.label)])
            
          } 
          
          #Randomly select one of the species
          spv2 <- sample(spv, 1)
          
          # Get the tip location for the sister sp.
          nodeX <- which(ctree$tip.label == spv2) 
          
         # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])

        }else if(ex$Type[j] == "RFG"){
          
          # Get all species within the family 
          fam <- dplyr::filter(jetz, BLFamilyLatin == ex$Sister_family[j])
          spv <- vector()
          
          for(x in 1:nrow(fam)){
            
            spv <- c(spv, paste0(fam$Genus[x], "_", fam$Species[x]) )
            
          } # Make the species group as a vector
          
          #Randomly select one of the species
          spv2 <- sample(spv, 1)
          
          # Get the tip location for the sister sp.
          nodeX <- which(ctree$tip.label == spv2) 
          
          # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])

        }else if(ex$Type[j] == "RCG"){
          
          # Get all species within the family 
          fam <- dplyr::filter(jetz, Clade == ex$Sister_clade[j])
          spv <- vector()
          
          for(x in 1:nrow(fam)){
            
            spv <- c(spv, paste0(fam$Genus[x], "_", fam$Species[x]))
            
          } # Make the species group as a vector
          
          #Randomly select one of the species
          spv2 <- sample(spv, 1)
          
          # Get the tip location for the sister sp.
          nodeX <- which(ctree$tip.label == spv2) 
          
          # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])

        }else if(ex$Type[j] == "ROG"){
          
          # Get all species within the family 
          ord <- dplyr::filter(jetz, order == ex$Sister_family[j])
          spv <- vector()
          
          for(x in 1:nrow(fam)){
            
            spv <- c(spv, paste0(ord$Genus[x], "_", ord$Species[x]) )
            
          } # Make the species group as a vector
          
          #Randomly select one of the species
          spv2 <- sample(spv, 1)
          
          # Get the tip location for the sister sp.
          nodeX <- which(ctree$tip.label == spv2) 
          
          # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])

        }
        
      } else if(ex$Type[j] %in% c("RSGG", "RSGG2")){
       ##' Scenario 3.1: Add species as a sister to a genus selected randomly from a supplied group of genera (RSGG) or a 
       ##'               random genus from a supplied family (RSGG2)
       
        if(ex$Type[j] == "RSGG"){
        
          #Split the supplied genera
          sp <- stringr::str_split(ex$Sister_genus[j], pattern = ";")
          
          #Randomly select one of the genera
          spv2 <- sample(sp[[1]], 1)
          
          if(length(ctree$tip.label[grep(paste0(spv2, "_"), ctree$tip.label)]) == 1){
            
            # Get the tip location for single species in the genus
            nodeX <- which(ctree$tip.label == ctree$tip.label[grep(paste0(spv2, "_"), ctree$tip.label)]) 
            
            # Bind the extinct sp.
            ctree <- bind(tree = ctree, node = nodeX,
                          per = PER, sp_name = ex$species[j])

          }else{
          
          # Get most recent common ancestor of genus 
          nodeX <- getMRCA(ctree, ctree$tip.label[grep(paste0(spv2, "_"), ctree$tip.label)]) 
          
          # Bind the extinct sp.
          ctree <- bind(tree = ctree, node = nodeX,
                        per = PER, sp_name = ex$species[j])

          }
          
        }
        
        if(ex$Type[j] == "RSGG2"){
          
          # Get all genera within the family 
          fam <- dplyr::filter(jetz, BLFamilyLatin == ex$Sister_family[j])
          spv <- vector()
          
          for(x in 1:nrow(fam)){
            
            spv <- c(spv, paste0(fam$Genus[x]) )
            
          } # Make the species group as a vector
          
          spv <- unique(spv)

          #Randomly select one of the genera
          spv2 <- sample(spv, 1)
          
          if(length(ctree$tip.label[grep(paste0(spv2, "_"), ctree$tip.label)]) == 1){
            
            # Get the tip location for single species in the genus
            nodeX <- which(ctree$tip.label == ctree$tip.label[grep(paste0(spv2, "_"), ctree$tip.label)]) 
         
            # Bind the extinct sp.
            ctree <- bind(tree = ctree, node = nodeX,
                          per = PER, sp_name = ex$species[j])

          } else{
            
            # Get most recent common ancestor of genus 
            nodeX <- getMRCA(ctree,ctree$tip.label[grep(paste0(spv2, "_"), ctree$tip.label)]) 
            
            # Bind the extinct sp.
            ctree <- bind(tree = ctree, node = nodeX,
                          per = PER, sp_name = ex$species[j])

          }
          
        } #eo if RSGG2

      } else {
        stop("Code / Type not recognised")
      }#eo main if statements

    } #eo for j
    
    return(ctree)          # Return the tree object
    
  }#eo for each
  
## Finish Tree ## - Should be 600 species added (10595) minus 13 already present (587 added, 10580 total)
class(ctreesComplete) <- "multiPhylo"    # Change the class
length(ctrees[[1]]$tip.label)            # Get the length of the tip labels (original)
length(ctreesComplete[[1]]$tip.label)    # Get the length of the tip labels (updated)
length(ctreesComplete[[1]]$tip.label) - length(ctrees[[1]]$tip.label) # Check the differences 

# Save the trees out
saveRDS(ctreesComplete, file = "JetzExtinct50_25082323.rds")
