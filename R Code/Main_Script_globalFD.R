
##ALL OPTIONS HAVE BEEN SET TO THOSE USED IN THE MAIN FD
#Analyses, but all can be changed.

##Use body size corrected traits (FALSE) or not (TRUE)
body_size_CORR <- FALSE

#Number of PCA axes to use in hypervolumes
no_axes <- 5

#body shape analysis (FALSE or TRUE); TRUE should only be used with 
#body size corrected = FALSE and no_axes = 4
body_shape <- FALSE

##FD or PD analysis
FD_PD <- "FD"

#trait_removal sensitivity analysis (NULL for FALSE, or name of trait to remove)
tra_remove <- NULL

#indiv imputation sensitivity analysis (NULL or replacement extinct2 object)
ind_impu <- NULL

#seabird sensitivity test (if TRUE, swap islandendemic column for 
#islandendemic2 column, which has the alternative seabird classifications)
seabird <- FALSE

##Source data and functions based on above selections
source("R Code\\Source_code_globalFD.R")


#######################################################################################
############################################
#########Main Analyses##################
#########################################################

##IF RUNNING ON A SINGLE MACHINE, USE the below code. However, the analyses for
#the paper were run on a HPC using array output, with the processing code outlined
#below this section.

cores = 8
cl = makeCluster(cores)
registerDoParallel(cl)
i = 1 #Dummy line for RStudio warnings

#####Hyper / convex / tree 
meth = "convex" #"hyper" or "convex" or "tree"
meth.hv = "svm" #"svm" or "gaussian"
gam = 0.5 #svm gamma parameter
GEOGz = "Yes" #"All" or "Yes" for island endemics
analyz = "main" #"main" or "beak"
dendroz = dendro_all#dendro_all or phylo_all[[1]]
use_phyloz = FALSE #TRUE for phylo, FALSE for dendrogram
rnz = 5#hypervolume replicate observed runs

#Number of future sims to use, set low here to run code
#more quickly.
#FD_all species = $Birdlife_all;
#FD_isl endemics = $Birdlife_isl
#PD_all species = $Birdtree_all;
#PD_isl endemics = $Birdtree_isl
#$Birdlife_isl2 and $Birdtree_isl2 = seabird sensitivity data
future_Fz = future$Birdlife_isl[1:2]#Future simulations (only using 2 for speed; 100 in main run)

##Observed values
main_run_all1_obs <- null_hyper(allSp2, 
                                future_F = future_Fz,
                                GEOG = GEOGz, analysis = analyz, 
                                axes_N = no_axes,
                                method = meth, observed = TRUE, 
                                method.hv = meth.hv, svmg = gam, 
                                BAND = BAND, SPP = SPP,  
                                dendro = dendroz, 
                                use_picante = TRUE,
                                use_phylo = use_phyloz,
                                rn = rnz)

#Percentage change in FD/PD across the four periods
perc_change(main_run_all1_obs[[1]]$Alpha)
perc_change(main_run_all1_obs[[1]]$Dispersion)#no dispersion for convex

#Observed Values
obs_vals <- main_run_all1_obs[[1]]

##Make the inset plots for these runs
if (meth == "hyper"){
  ggtitlz <- paste0(GEOGz, " - ", meth.hv, gam, "_BS", body_size_CORR) 
} else{
  ggtitlz <- paste0(GEOGz, " - ", meth, "_BS", body_size_CORR)
}

IP1 <- hyper_plot(obs_vals, ggtitl = ggtitlz , method = meth,
                  inset_only = TRUE)

jpeg(file = paste0(ggtitlz,".jpeg"), width = 20, height = 10, 
     units = "cm", res = 300)
IP1
dev.off()


##Null model values

n1 = 10 #number of null iterations (Set low here to run code fast)

main_run_all1_null = foreach(i=seq(from=1, to=n1, by=1), 
                             .inorder = FALSE)  %dopar% {
                               library(dplyr)
                               library(BAT)
                               #observed needs to = FALSE for null run
                               null_hyper(allSp2, GEOG = GEOGz, analysis = analyz,
                                          future_F = future_Fz,
                                          axes_N = no_axes, method = meth, 
                                          observed = FALSE,
                                          method.hv = meth.hv, svmg = gam,
                                          BAND = BAND, SPP = SPP,
                                          dendro = dendroz,
                                          use_picante = TRUE,
                                          use_phylo = use_phyloz)
                             }
#Get ES values and P-values
main_run_all1_ZP <- null_hyper_ZP(obs_vals, 
                                  main_run_all1_null,
                                  plot_null = FALSE)

main_run_all1 <- list(obs_vals, main_run_all1_null, 
                      main_run_all1_ZP)
#save(main_run_all1, file = "main_run_all1.R")

##################################################
######PROCESSING THE ARRAY OUTPUT#################
##################################################

###For the run analyses, array processing was used on a HPC.

#The output files are saved in the Results directory, with
#sub-directories for each separate analysis:
#All = all species; Isl = island endemics
#0.5 = hypervolumes; convex = convex hulls; tree = dendrogram
#phylo = phylogeny
#Isl2 = alternative island endemic analysis run

##in each output file, the 1st element of the list is
#a list with two elements: (1) the observed values for all time periods,
#including the median values for the future period, 
#and (2) the individual values for each of the N future runs. 
#The subsequent elements of the list are then all 
#individual null model runs.

gw <- getwd()
setwd("Results\\main results")
DL <- list.dirs(getwd(), recursive = FALSE)

main_res_list <- vector("list", length = length(DL))
phylo_res_list <- vector("list", length = 2)
k <- 1

for (j in 1:length(DL)){
  
  wd <- DL[j]
  
  main_run_1 <- vector("list", length = 4)
  
  setwd(wd)
  lf_isl <- list.files(wd)
  
  main_run_0 <- readRDS(lf_isl)
  
  #for the phylo results, we have a list with elements for each of
  #the phylogenies, so have to extract each with lapply. But to do
  #this, we store it in a new list, and then delete all but the first
  #tree's results, so that main_run_1 only has the results for one tree
  #to enable it to work with the plotting function
  if (grepl("phylo", wd)){
    
    main_run_2 <- vector("list", length = 4)
    
    #observed values
    main_run_2[[1]] <-  lapply(main_run_0, 
                               function(x) x[[1]][[1]])
    
    #individual future runs values
    main_run_2[[4]] <-  lapply(main_run_0, 
                               function(x) x[[1]][[2]])
    
    #null
    main_run_2[[2]] <- lapply(main_run_0, 
                              function(x) x[2:length(x)])
    
    #remember, null_hyper_ZP changes the observed values (main_run_1) to
    #be H,C,C instead of P,H,C
    main_run_2[[3]] <- lapply(main_run_0, function(x){
      null_hyper_ZP(x[[1]][[1]], x[2:length(x)], plot_null = FALSE)})
    
    phylo_res_list[[k]] <- main_run_2
    names(phylo_res_list)[k] <- sub(".*ArrayOut", "", wd)
    k <- k + 1
    
    #just store the random tree for processing below
    if ((grepl("Isl2", wd))){
      main_run_0 <- main_run_0[[20]]
    } else {
      main_run_0 <- main_run_0[[31]]
    }
    
    
  }#eo if phylo
  
  #observed values
  main_run_1[[1]] <-  main_run_0[[1]][[1]]
  #individual future runs values
  main_run_1[[4]] <-  main_run_0[[1]][[2]]
  
  #null
  main_run_1[[2]] <- main_run_0[2:length(main_run_0)]
  
  #remember, null_hyper_ZP changes the observed values (main_run_1) to
  #be H,C,C instead of P,H,C
  main_run_1[[3]] <- null_hyper_ZP(main_run_1[[1]], main_run_1[[2]],
                                   plot_null = FALSE)
  
  #save with specific name for this analysis
  main_res_list[[j]] <- main_run_1
  names(main_res_list)[j] <- sub(".*ArrayOut", "", wd)
  
}#eo j

setwd(gw)

#save(main_res_list, file = "main_res_list.R")
#save(phylo_res_list, file = "phylo_res_list.R")

###########################################################
######Plotting Main Null Model Analyses############
##########################################################

##build multi plots
AJ <- list(c("All_0.5", "a) All - hypervolume0.5", "hyper", "nrd"), 
           c("Isl_0.5", "b) Island - hypervolume0.5", "hyper", "nrd"),
           c("All_convex", "a) All - FRic", "convex", "nrd"), 
           c("Isl_convex", "b) Island - FRic", "convex", "nrd"),
           c("All_tree", "a) All - FD", "tree", "nrd"), 
           c("Isl_tree", "b) Island - FD", "tree", "nrd"),
           c("All_phylo", "a) All - phylo", "tree", "nrd"),
           c("Isl_phylo", "b) Island - phylo", "tree", "nrd"),
           c("Isl2_0.5", "b) Island2 - hypervolume0.5", "hyper", "nrd"),
           c("Isl2_phylo", "b) Island2 - phylo", "tree", "nrd"))

main_plots_comb <- lapply(AJ, 
                          function(x) hyper_plot(beak_run_F = main_res_list[[x[1]]],
                                                 ggtitl = x[2], method = x[3],
                                                 side = "one", BW = x[4]))

#inset plots (plot all together)
all_inset <- lapply(main_plots_comb, function(x) x[[2]])

jpeg(file = "main_inset.jpeg", width = 10, height = 43, 
     units = "cm", res = 300)
gridExtra::grid.arrange(grobs = all_inset, nrow = 10)
dev.off()

##HYPERVOLUME 0.5###
hyper0.5_plots <- lapply(main_plots_comb, function(x) x[[1]])[1:2]

jpeg(file = "hyper0.5.jpeg", width = 16, height = 14, units = "cm", res = 300)
gridExtra::grid.arrange(grobs = hyper0.5_plots, nrow = 2)
dev.off()

##CONVEX - FRIC###
convex_plots <- lapply(main_plots_comb, function(x) x[[1]])[3:4]

jpeg(file = "convex.jpeg", width = 16, height = 14, units = "cm", res = 300)
gridExtra::grid.arrange(grobs = convex_plots, nrow = 2)
dev.off()

##TREE - FD###
tree_plots <- lapply(main_plots_comb, function(x) x[[1]])[5:6]

jpeg(file = "tree.jpeg", width = 16, height = 14, units = "cm", res = 300)
gridExtra::grid.arrange(grobs = tree_plots, nrow = 2)
dev.off()

##TREE - PD###
phylo_plots <- lapply(main_plots_comb, function(x) x[[1]])[7:8]

jpeg(file = "phylo.jpeg", width = 16, height = 14, units = "cm", res = 300)
gridExtra::grid.arrange(grobs = phylo_plots, nrow = 2)
dev.off()

##Seabird sensitivity - hyper and phylo###
#NOTE - for seabird > island > alpha > All->Cur, the diamond
#colour needs to be manually changed to red as it is based on 
#majority rules, and only 1/20 phylogenies is significant for this;
#it just happens to be the one used for plotting.
seabird_plots <- lapply(main_plots_comb, function(x) x[[1]])[9:10]

jpeg(file = "seabird.jpeg", width = 16, height = 14, units = "cm", res = 300)
gridExtra::grid.arrange(grobs = seabird_plots, nrow = 2)
dev.off()


##################################################
###Results table######
################################################################

##FD

#need to remove phylo and analyse them separately
main_res_list2 <- main_res_list[c("All_0.5", "All_convex", 
                                  "All_tree",
                                  "Isl_0.5","Isl_convex",
                                  "Isl_tree",
                                  "Isl2_0.5")]

res_table <- matrix(NA, nrow = length(main_res_list2), 
                    ncol = 40)
rownames(res_table) <- c("All_0.5", 
                         "Isl_0.5",
                         "All_convex", 
                         "Isl_convex",
                         "All_tree", 
                         "Isl_tree",
                         "Isl2_0.5")
colnames(res_table) <- c("P_Alpha", "H_Alpha","C_Alpha",
                         "F_Alpha",
                         "PH_Alpha_Perc",
                         "PC_Alpha_Perc",
                         "HC_Alpha_Perc",
                         "CF_Alpha_Perc",
                         "PC_Alpha_Z","HC_Alpha_Z","CF_Alpha_Z", 
                         "PC_Alpha_P","HC_Alpha_P",
                         "CF_Alpha_P",
                         "F_Alpha_RL", "F_Alpha_RU", "F_Alpha_Med", 
                         "F_Alpha_Perc_RL", "F_Alpha_Perc_RU", 
                         "F_Alpha_Perc_Med", 
                         "P_Disp", "H_Disp","C_Disp",
                         "F_Disp",
                         "PH_Disp_Perc",
                         "PC_Disp_Perc","HC_Disp_Perc",
                         "CF_Disp_Perc",
                         "PC_Disp_Z",
                         "HC_Disp_Z","CF_Disp_Z",
                         "PC_Disp_P",
                         "HC_Disp_P","CF_Disp_P",
                         "F_Disp_RL", "F_Disp_RU", "F_Disp_Med", 
                         "F_Disp_Perc_RL", "F_Disp_Perc_RU", 
                         "F_Disp_Perc_Med")


for (i in 1:length(rownames(res_table))){
  
  k <- main_res_list2[[rownames(res_table)[i]]]
  
  res_table[i, 1:4] <- round(k[[1]]$Alpha,2)
  res_table[i, 5:8] <- round(perc_change(k[[1]]$Alpha),2)
  res_table[i, 9:11] <- filter(k[[3]], Metric == "Alpha")$Z
  res_table[i, 12:14] <- filter(k[[3]], Metric == "Alpha")$P
  if (!identical(filter(k[[3]], Metric == "Alpha")$Period, 
                 c("All", "IUCN", "Fut"))){
    stop ("periods out of order")
  }
  ##future sim runs
  res_F_alp <- unlist(k[[4]][["alpha"]])
  res_table[i, 15:16] <- round(range(res_F_alp),2)
  res_table[i, 17] <- round(median(res_F_alp),2)
  if (!round(median(res_F_alp),3) == 
      round(k[[1]]$Alpha[4],3)) stop("fAINT")
  #work out percentage change between observed C value,
  #and each of the F values. Only the first value out of
  #perc_change is taken, others ignored
  res_F_alp_per <- sapply(res_F_alp, 
                          function(x) perc_change(c(k[[1]]$Alpha[3], x, 
                                                    99,99)))[1,]
  res_table[i, 18:19] <- round(range(res_F_alp_per),2)
  #median of perc change values = the percentage change calculated
  #using the median F_alpha value.
  res_table[i, 20] <- round(median(res_F_alp_per),2)
  
  ##dont need to fill in the dispersion metric for convex as only alpha
  if (!rownames(res_table)[i] %in% c("All_convex", "Isl_convex")){
    CNR <- "Dispersion" 
    
    res_table[i, 21:24] <- round(k[[1]][,CNR],2)
    res_table[i, 25:28] <- round(perc_change(k[[1]][,CNR]),2)
    res_table[i, 29:31] <- filter(k[[3]], Metric == CNR)$Z
    res_table[i, 32:34] <- filter(k[[3]], Metric == CNR)$P
    if (!identical(filter(k[[3]], Metric == CNR)$Period, 
                   c("All", "IUCN", "Fut"))){
      stop ("periods out of order")
    }
    ##future sim runs
    res_F_dis <- unlist(k[[4]][["dispersion"]])
    res_table[i, 35:36] <- round(range(res_F_dis),2)
    res_table[i, 37] <- round(median(res_F_dis),2)
    if (!round(median(res_F_dis),3) == 
        round(k[[1]]$Dispersion[4],3)) stop("fAINT2")
    res_F_dis_per <- sapply(res_F_dis, 
                            function(x) perc_change(c(k[[1]]$Dispersion[3], x, 
                                                      99,99)))[1,]
    res_table[i, 38:39] <- round(range(res_F_dis_per),2)
    res_table[i, 40] <- round(median(res_F_dis_per),2)
  }#eo if convex
  
}#eo for i

res_table <- as.data.frame(res_table)

res_table <- apply(res_table, 2, as.numeric)
rownames(res_table) <- c("All_0.5", 
                         "Isl_0.5",
                         "All_convex", 
                         "Isl_convex",
                         "All_tree", 
                         "Isl_tree",
                         "Isl2_0.5")

##Make columns for the future and C2F median values, 
#with the range in parentheses
sasg1 <- apply(res_table, 1, function(x){
  paste0(x["F_Alpha_Med"],
         " (", round(x["F_Alpha_RL"],1),
         ":", round(x["F_Alpha_RU"],1), ")")
})
sasg2 <- apply(res_table, 1, function(x){
  paste0(x["F_Alpha_Perc_Med"],
         " (", round(x["F_Alpha_Perc_RL"],1),
         ":", round(x["F_Alpha_Perc_RU"],1), ")")
})
write.csv(sasg1)#F Alpha values
write.csv(sasg2)#C2F perc change values

write.csv(res_table, file = "main_res_table.csv")

##PD
prl_All <- PD_multi_format(phylo_res_list$All_phylo)#all
prl_isl <- PD_multi_format(phylo_res_list$Isl_phylo)#isl
prl_isl2 <- PD_multi_format(phylo_res_list$Isl2_phylo)#isl2

if (!identical(names(prl_All[[1]]), 
               names(prl_isl[[1]]))) stop("MOS")

if (!identical(names(prl_All[[1]]), 
               names(prl_isl2[[1]]))) stop("MOS")

pd_res_table <- bind_rows(prl_All[[1]], prl_isl[[1]],
                          prl_isl2[[1]]) %>%
  as.data.frame() %>%
  round(2)
rownames(pd_res_table) <- c("All_phylo", "Isl_phylo",
                            "Isl2_phylo")

write.csv(pd_res_table, file = "pd_res_table.csv")


##################################################
###Temporal trend plot: using percentage change######
################################################################

#NB: Source code has to be run with FD_PD == "FD"

res_table <- as.data.frame(res_table)

#Remove seabird sensitivity results from PD results table
pd_res_table <- pd_res_table[which(rownames(pd_res_table) %in% 
                                     c("All_phylo", "Isl_phylo")),]

###Percentage changes are all relative to P
#i.e. P2H, P2C, P2F

##subset main results tables
ttp_df_FD <- res_table %>%
  select("PH_Alpha_Perc", "PC_Alpha_Perc",
         "P_Alpha", "F_Alpha") %>%
  filter(rownames(.) %in% c("All_0.5", "Isl_0.5")) %>%
  mutate(across(1:4, as.numeric))

ttp_df2_FD <- ttp_df_FD %>%
  rowwise() %>%
  mutate("PF_Alpha_Perc" = 
           perc_change(c(P_Alpha, F_Alpha, 99, 99))[1]) %>%
  select("PH_Alpha_Perc", "PC_Alpha_Perc", 
         "PF_Alpha_Perc") %>%
  mutate(across(3, as.numeric)) %>% 
  as.data.frame() 
rownames(ttp_df2_FD) <- rownames(ttp_df_FD)

ttp_df_PD <- pd_res_table %>%
  select("PH_Alpha_Perc" = "PH_Alpha_Perc_Med", 
         "PC_Alpha_Perc" = "PC_Alpha_Perc_Med",
         "PF_Alpha_Perc" = "PF_Alpha_Perc_Med") %>%
  mutate(across(1:3, as.numeric))


ttp_df <- rbind.data.frame(ttp_df2_FD, ttp_df_PD)
ttp_df2 <- as.data.frame(ttp_df)
rownames(ttp_df2) <- vapply(rownames(ttp_df), 
                            function(x){
                              switch(x, 
                                     "All_0.5" = "FD_all",
                                     "Isl_0.5" = "FD_isl",
                                     "All_phylo" = "PD_all",
                                     "Isl_phylo" = "PD_isl")
                            }, FUN.VALUE = character(1))


##convert results data in figure data:
#P = 100%, with H, C and F % decreases from
#this.
#y-axis = % persisting FD
ttp_df3 <- matrix(nrow = 4, ncol = 4)
colnames(ttp_df3) <- c("P", "H", "C", "F")
rownames(ttp_df3) <- rownames(ttp_df2)
ttp_df3[,1] <- 100

for (i in 1:nrow(ttp_df3)){
  #we use + as the perc change values are negative;
  #and all are relative to 100%
  ttp_df3[i,2] <- 100 + 
    ttp_df2[i,"PH_Alpha_Perc"]
  
  ttp_df3[i,3] <- 100 + 
    ttp_df2[i,"PC_Alpha_Perc"]
  
  ttp_df3[i,4] <-  100 +  
    ttp_df2[i,"PF_Alpha_Perc"]
}

ttp_df3 <- as.data.frame(ttp_df3) %>%
  mutate("Metric" = rownames(ttp_df3))

##add in species richness values (note, this needs to
#be run with the BirdLife data)
#ALL SPECIES
AS1P <- nrow(allSp2)
AS1H <- nrow(allSp2) - table(allSp2$ExtinctionPeriod)["All"]
AS1C <- filter(allSp2, Extant == "Extant") %>% nrow()
AS1F <- AS1C - length(future[[1]][[1]])
AS1 <- c(AS1P, AS1H, AS1C, AS1F) %>%
  as.vector()
AS1P1 <- perc_change(AS1)
AS1P2 <- perc_change(c(AS1P, AS1F, 999,999))
ASR1 <- data.frame("P" = 100, "H" = (100 + AS1P1[1]),
                   "C" = (100 + AS1P1[2]), 
                   "F" = (100 + AS1P2[1]), "Metric" = "SR_all")
rownames(ASR1) <- c("SR_all")

#ISL ENDEMICS
ASallSp2 <- filter(allSp2, IslandEndemic == "Yes")
AS2P <- nrow(ASallSp2)
AS2H <- nrow(ASallSp2) - table(ASallSp2$ExtinctionPeriod)["All"]
AS2C <- filter(ASallSp2, Extant == "Extant") %>% nrow()
AS2F <- AS2C - length(future[[2]][[1]])
AS2 <- c(AS2P, AS2H, AS2C, AS2F) %>%
  as.vector()
AS2P1 <- perc_change(AS2)
AS2P2 <- perc_change(c(AS2P, AS2F, 999,999))
ASR2 <- data.frame("P" = 100, "H" = (100 + AS2P1[1]),
                   "C" = (100 + AS2P1[2]), 
                   "F" = (100 + AS2P2[1]), "Metric" = "SR_isl")
rownames(ASR2) <- c("SR_isl")

ttp_df3b <- rbind.data.frame(ttp_df3, ASR1, ASR2)

ttp_df4 <- tidyr::pivot_longer(data = ttp_df3b, 
                               cols = P:'F')
ttp_df4$name <- factor(ttp_df4$name, 
                       levels = c("P", "H", "C", "F"))

gg_ttp <- ggplot(data = ttp_df4, aes(name, value)) +
  geom_point(aes(col = Metric, shape = Metric), size = 3) +
  scale_color_manual(values =
                       c("#48ABAE", "#48ABAE",
                         "#F3D54D", "#F3D54D",
                         "Black", "Black")) +
  scale_shape_manual(values = c(17, 16, 17,
                                16, 17, 16)) +

  ylab("Percentage of diversity remaining") +
  xlab("Time period") + 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  labs(col = NULL, shape = NULL) +
  geom_line(data=filter(ttp_df4, 
                        name %in% c("P", "H", "C")),
            aes(col = Metric, group = Metric)) +
  geom_line(data=filter(ttp_df4, 
                        name %in% c("C","F")),
            linetype = 2, aes(col = Metric, group = Metric)) +
  scale_x_discrete(labels = c(expression('130K BP'["[All]"]), 
                              expression('1500 CE'["[IUCN]"]), 
                              "Current", "Future"))


jpeg(file = "temporal_trend.jpeg", width = 12, 
     height = 8, units = "cm", res = 300)
gg_ttp
dev.off()

#######################################################
###########INDIVIDUAL TRAIT ANALYSES##########################
#######################################################
  
#NULL ITERATIONS SET TO 99 here for SPEED

#Note - SES p-values presented here are two-tailed;
#they need to be converted to one-tailed for comparison
#with values from the ES approach

##Set up analyses for all species and island endemics,
#focused on median and SD trait values
GT <- list(c("All", "a) All", "Median", 0.1, 0.05,
             names(future)[1]),
           c("Yes", "b) Island", "Median", 1, 0.1,
             names(future)[2]),
           c("All", "c) All", "SD", 500, 0.01,
             names(future)[1]),
           c("Yes", "d) Island", "SD", 500, 0.1,
             names(future)[2]))

###for mass 
Mass_plots_comb <- lapply(GT, function(x) null_plot(future_F = future[[x[6]]],
                                                    TRAIT = "Mass", 
                                                    GEOG = x[1],
                                               TEST = x[3], 
                                               n = 99,
                                 titl = x[2], 
                                 bandwidth = as.numeric(x[4]),
                               #  color = x[6], 
                                 hist = FALSE))

Mass_plots <- lapply(Mass_plots_comb, function(x) x[[1]][[1]])
Mass_inset <- lapply(Mass_plots_comb, function(x) x[[1]][[2]])

jpeg(file = "Mass.jpeg", width = 25, height = 23, units = "cm", res = 300)
gridExtra::grid.arrange(grobs = Mass_plots, nrow = 2)
dev.off()

jpeg(file = "Mass_inset.jpeg", width = 10.5, height = 7, units = "cm", res = 300)
gridExtra::grid.arrange(grobs = Mass_inset, nrow = 2)
dev.off()

#% change between the four periods
Mass_perc_change <- lapply(Mass_plots_comb, function(x) x[[2]])[1:2]#as 3 and 4 just repeat
Mass_PC <- lapply(Mass_perc_change, function(x){
  #v2 new, v1 original  
  perc <- function(v1, v2){((v2 - v1) / v1) * 100}   
  x2 = x[,1:2]
  apply(x2, 2, function(y){
    p2 <- perc(y[1], y[3]) #P2C
    p3 <- perc(y[2], y[3]) #H2C
    p4 <- perc(y[3], y[4]) #C2F
    c("All" = p2, "IUCN" = p3, "fut" = p4)
  })
})
Mass_PC

##Observed values
lapply(Mass_plots_comb, function(x) x[[2]])
## Null model results (all median, isl median, all sd, isl sd)
Mass_ES <- lapply(Mass_plots_comb, function(x) x[[3]])
Mass_ES

##Future values (from individual runs); remember, for
#main results and plots, it uses the median of these
lapply(Mass_plots_comb, function(x) x[[4]])
#range & median together
lapply(Mass_plots_comb, function(x){
  frf <- apply(x[[4]], 1, range)
  frf2 <- apply(x[[4]], 1, median)
  rbind(frf, "Median" = frf2)
})[1:2] #as 3 and 4 repeat

###for hWI 
HWI_plots_comb <- lapply(GT[1:2], 
                         function(x) null_plot(future_F = future[[x[6]]],
                                               TRAIT = "Hand.Wing.Index", 
                                               GEOG = x[1],
                                               TEST = x[3], 
                                               n = 99,
                                               titl = x[2], 
                                               bandwidth = as.numeric(x[5]),
                                             #  color = x[6], 
                                              hist = FALSE))

HWI_plots <- lapply(HWI_plots_comb, function(x) x[[1]][[1]])
HWI_inset <- lapply(HWI_plots_comb, function(x) x[[1]][[2]])


jpeg(file = "HWI.jpeg", width = 25, height = 23, units = "cm", res = 300)
gridExtra::grid.arrange(grobs = HWI_plots, nrow = 2)
dev.off()

jpeg(file = "HWI_inset.jpeg", width = 10.5, height = 7, units = "cm", res = 300)
gridExtra::grid.arrange(grobs = HWI_inset, nrow = 2)
dev.off()

#% change between the four periods
HWI_perc_change <- lapply(HWI_plots_comb, function(x) x[[2]])[1:2]
HWI_PC <- lapply(HWI_perc_change, function(x){
  #v2 new, v1 original  
  perc <- function(v1, v2){((v2 - v1) / v1) * 100}   
  x2 = x[,1:2]
  apply(x2, 2, function(y){
    p2 <- perc(y[1], y[3]) #P2C
    p3 <- perc(y[2], y[3]) #H2C
    p4 <- perc(y[3], y[4]) #C2F
    c("All" = p2, "IUCN" = p3, "Fut" = p4)
  })
})
HWI_PC

##Observed values
lapply(HWI_plots_comb, function(x) x[[2]])
## Null model results
HWI_ES <- lapply(HWI_plots_comb, function(x) x[[3]])
HWI_ES

##Future values (from individual runs); remember, for
#main results and plots, it uses the median of these
lapply(HWI_plots_comb, function(x) x[[4]])
#range & median together
lapply(HWI_plots_comb, function(x){
  frf3 <- apply(x[[4]], 1, range)
  frf4 <- apply(x[[4]], 1, median)
  rbind(frf3, "Median" = frf4)
})[1:2] #as 3 and 4 repeat



 #########Beak Analyses##################

 cores = 8
 cl = makeCluster(cores)
 registerDoParallel(cl)
 i = 1 #Dummy line for RStudio warnings
 
 
 #Null iterations set to 9 HERE FOR SPEED
 n2 = 9
 gam = 0.5
 meth.hv = "svm" #svm or gaussian

 beak_res <- vector("list", length = 2)
 k <- 1
 
 #GEOG = "All", "Yes" (for island endemics) or "No" (For continental)
 for (j in c("All", "Yes")){
   
   #Number set low for SPEED
   if (j == "All"){
     future_Fz <- future$Birdlife_all[1:2]
   } else {
     future_Fz <- future$Birdlife_isl[1:2]
   }
   
   #####Hyper
   #rn set to 5 for SPEED
   beak_res[[k]][[1]] <- null_hyper(allSp2, 
                                    future_F = future_Fz,
                                    GEOG = j, 
                                    analysis = "beak", 
                                    method = "hyper", 
                                    method.hv = meth.hv,
                                    observed = TRUE, svmg = gam,
                                    BAND = BAND, SPP = SPP,
                                    dendro = dendro_all,
                                    rn = 2)
   
   beak_res[[k]][[2]] = foreach(i=seq(from=1, to=n2, by=1), 
                                .inorder = FALSE)  %dopar% { 
     library(dplyr)
     library(BAT)
     #observed needs to = FALSE for null run
     null_hyper(allSp2, future_F = future_Fz,
                GEOG = j, analysis = "beak", 
                method = "hyper", observed = FALSE, 
                method.hv = meth.hv, svmg = gam,
                BAND = BAND, SPP = SPP,
                dendro = dendro_all)
   }
   
   k <- k + 1
 }#eo j
 
 names(beak_res) <- c("All", "Isl")
 
 beak_res_0.5 <- beak_res

 beak_resP <- beak_res_0.5
 
 ##plots
 beak_ti <- c("All", "Isl")
 
 ##Change hist to TRUE to return histograms of
 #null model values for visual checks
 beak_plots_comb <- lapply(beak_ti, function(x) {
   
   beak_ti2 <- switch(x,
                       "All" = "a) All",
                       "Isl" = "b) Island")
   
   null_plot(titl = paste0(beak_ti2,"_0.5"), 
        #     color =  "#889eb9",
             analysis = "beak", 
             beak_data = beak_resP[[x]],
             hist = FALSE)
   
 })
 
 beak_plots <- lapply(beak_plots_comb, function(x) x[[1]][[1]])
 beak_inset <- lapply(beak_plots_comb, function(x) x[[1]][[2]])
 
 
 jpeg(file = "beak_0.5.jpeg", width = 25, height = 10, units = "cm", res = 300)
 gridExtra::grid.arrange(grobs = beak_plots, nrow = 1)
 dev.off()
 
 jpeg(file = "beak_inset_0.5.jpeg", width = 7, height = 8, units = "cm", res = 300)
 gridExtra::grid.arrange(grobs = beak_inset, nrow = 2)
 dev.off()
 
 
 #% change between the four periods
 beak_perc_change <- lapply(beak_plots_comb, function(x) x[[2]])[1:2]
 beak_PC <- lapply(beak_perc_change, function(x){
   #v2 new, v1 original  
   perc <- function(v1, v2){((v2 - v1) / v1) * 100}   
   y = x[,1]
   p2 <- perc(y[1], y[3]) #P2C
   p3 <- perc(y[2], y[3]) #H2C
   p4 <- perc(y[3], y[4]) #C2F
   c("All" = p2, "IUCN" = p3, "Fut" = p4)
 })
 beak_PC

 ##Observed values
 lapply(beak_plots_comb, function(x) x[[2]])
 ## Null model results
 beak_ES <- lapply(beak_plots_comb, function(x) x[[3]])
 beak_ES

###combine all single trait results in big table for SI
 #Percentage change
 SI_PC <- matrix(nrow = 6, ncol = 4)
 rownames(SI_PC) <- rep(rownames(Mass_PC[[1]]), 2)
 colnames(SI_PC) <- c("Mass_median", "Mass_SD", 
                      "HWI_median", "Beak")
 SI_PC[,1:2] <- do.call(rbind, Mass_PC)
 SI_PC[,3] <- do.call(rbind, HWI_PC)[,1]
 SI_PC[,4] <- do.call(c, beak_PC)
 round(SI_PC, 2) %>% write.csv()
 
 #ES / SES values
 SI_ES_mass <- matrix(nrow = 6, ncol = 16)
 rownames(SI_ES_mass) <- rep(rownames(Mass_PC[[1]]), 2)
 colnames(SI_ES_mass) <- c("Mass_Median_ES", "Mass_Median_ES_P",
                           "Mass_Median_SES", "Mass_Median_SES_P",
                           "Mass_SD_ES", "Mass_SD_ES_P",
                           "Mass_SD_SES", "Mass_SD_SES_P",
                           "HWI_ES", "HWI_ES_P",
                           "HWI_SES", "HWI_SES_P",
                           "beak_ES", "beak_ES_P",
                           "beak_SES", "beak_SES_P")
 SI_ES_mass[1:3,1:4] <- Mass_ES[[1]]
 SI_ES_mass[4:6,1:4] <- Mass_ES[[2]]
 SI_ES_mass[1:3,5:8] <- Mass_ES[[3]]
 SI_ES_mass[4:6,5:8] <- Mass_ES[[4]]
 SI_ES_mass[1:3,9:12] <- HWI_ES[[1]]
 SI_ES_mass[4:6,9:12] <- HWI_ES[[2]]
 SI_ES_mass[1:3,13:16] <- beak_ES[[1]]
 SI_ES_mass[4:6,13:16] <- beak_ES[[2]]
 
 t(SI_ES_mass)
 
 
 ##################################################################
 ###########CONTRIBUTION ANALYSIS######################
 #############################################################
 
 #NOTE - this takes ~1 hour to run for each data set (e.g. all species, FD),
 #but the results can be loaded in separately (see below)
 
 cores = 8
 cl = makeCluster(cores)
 registerDoParallel(cl)
 i = 1 #Dummy line for RStudio warnings
 #stopCluster(cl)
 
#All species or island endemic
 allSp2_contr <- allSp2#'allSp2' (all species) or 'filter(allSp2, IslandEndemic == "Yes")'

 #FD_PD = "FD" or "PD" for dendrogram/trees
 #hyp_meth = "svm" or "gaussian"
 #type = "ES" or "SES"
 #rel = return relative contribution values or not
 #n_contr = number of replications for FD_hyper
 #null_contr_N = number of null iterations
 
 n1 = 1 #1 for FD and 50 for PD
 
 tree_contr <- list(dendro_all) #'list(dendro_all)' or 'phylo_all'
 
 k5Cont = foreach(i=seq(from=1, to=n1, by=1))  %dopar% {
   library(dplyr)
   library(BAT)
   library(ape)
   contribution_own(allSp_data = allSp2_contr, 
                    hyp_meth = "svm",
                    SPP = SPP,
                    rel = FALSE, axes_N = 5, 
                    FD_PD = FD_PD, 
                    tree = tree_contr[[i]], 
                    type = "ES",
                    null_contr_N = 999)
   
 }
 
 ######Final plot for paper################
 
 #Results of the contribution analyses are saved in the Results 
 #directory and are loaded below

 cont_p4_FDPD <- vector("list", length = 2)
 
 #1 = FD, 2 = PD
 for (i in 1:2){
   
   if (i == 1){
     CFDPD <- "FD"
   } else if (i == 2) {
     CFDPD <- "PD"
   }
   
   #All species
   if (i == 1){
     load("Results\\Contribution\\contributions_relF_FD_all.R")
   } else if (i == 2){
     load("Results\\Contribution\\contributions_relF_PD_all.R")
     
   }
   cont_p4_all <- lapply(k5Cont, function(x) x[[1]])
   #Randomly using 1st set of contribution values
   cont_p4_all <- cont_p4_all[[1]][[1]]
   cont_p4_all$Contribution <- c(rep(paste0(CFDPD,"_all"), 5), 
                                 rep("SR_all", 5))
   #Island endemics
   if (i == 1){
     load("Results\\Contribution\\contributions_relF_FD_isl.R")
   } else if (i == 2){
     load("Results\\Contribution\\contributions_relF_PD_isl.R")
   }
   cont_p4_isl <- lapply(k5Cont, function(x) x[[1]])
   cont_p4_isl <- cont_p4_isl[[1]][[1]]
   cont_p4_isl$Contribution <- c(rep(paste0(CFDPD,"_isl"), 5), 
                                 rep("SR_isl", 5))
   
   cont_p4 <- bind_rows(cont_p4_all, cont_p4_isl)
   cont_p4$Contribution <- factor(cont_p4$Contribution, 
                                  levels = c(paste0(CFDPD,"_all"), 
                                             "SR_all",
                                             paste0(CFDPD,"_isl"), 
                                             "SR_isl"))
   
   #Remove SR now,
   cont_p4 <- filter(cont_p4, 
                     Contribution %in% c("FD_all", "FD_isl",
                                         "PD_all", "PD_isl"))
   
   cont_p4_FDPD[[i]] <- cont_p4
   
 }#eo for i
 
 
 cont_p4_FDPD2 <- as.data.frame(do.call(rbind, cont_p4_FDPD))
 
 cont_p4_FDPD2$ExtinctionType <- factor(cont_p4_FDPD2$ExtinctionType,
                                        levels = c("NT", "TH", 
                                                   "IUCN", "EPU", "EPA",
                                                   "LR"))
 cont_p4_FDPD2$ExtinctionType[which(cont_p4_FDPD2$ExtinctionType == "NT")] <- "LR"
 
 cont_p4_FDPD2$ExtinctionType <- factor(cont_p4_FDPD2$ExtinctionType,
                                        levels = c("LR", "TH", "IUCN", "EPA", "EPU"))
 
 
 
 gg_con <- ggplot(data = cont_p4_FDPD2) + 
   geom_bar(position="fill", stat="identity",
            aes(y = n, x = Contribution, 
                fill = ExtinctionType)) +
   theme(panel.background = element_blank(), 
         axis.line = element_line(colour = "black")) +
   theme(axis.text = element_text(size=10),
         axis.title = element_text(size = 13),
         legend.text=element_text(size=11)) +
   scale_fill_manual(values = c("#BDC2E3","#4A6FE3", "#E6B8C0", 
                                "#E18095", "#7F000D"),
                     labels = c("LR", "TH", 
                                expression(EX["IUCN"]),
                                expression(EX["A"]), 
                                expression(EX["U"]))) +
   theme(legend.text.align = 0) +
   xlab("Contribution type") + 
   ylab("Proportion of total diversity") +
   labs(fill = "")
 
 jpeg(filename = "Figure4_contr.jpeg", 
      width = 15, height = 12, units = "cm", res = 300)
 gg_con
 dev.off()
 

#############################################################
########Morphospace Figure 1###############################
###########################################################
morpho_plot <- allSp2 %>%
  select(species, PC1, PC2, status)

#Stored FD contribution values (dendrogram) for all species
load("Results\\Contribution\\contributions_relF_FD_all.R")

contrVals_p4 <- lapply(k5Cont, function(x) x[[3]])

cdf <- data.frame("species" = names(contrVals_p4[[1]]),
                  "contribution" = as.vector(contrVals_p4[[1]]))

cdf <- cdf[match(morpho_plot$species,
                  cdf$species),]

if (!identical(morpho_plot$species, cdf$species)){
  stop("the Band")
}

morpho_plot$contribution <- cdf$contribution

IUCN <- morpho_plot$status ; unique(IUCN)
IUCN <- ifelse(IUCN == "LC" | IUCN == "NT" | IUCN == "DD", "LR", IUCN)
IUCN <- ifelse(IUCN == "CR" | IUCN == "EN" | IUCN == "VU", "TH", IUCN)

morpho_plot$IUCN <- IUCN

morpho_plot$IUCN <- factor(morpho_plot$IUCN, 
                           levels=c("LR", "TH", "EX"))
morpho_plot$contributionPer <- (morpho_plot$contribution/sum(morpho_plot$contribution)) * 100

g1 <- ggplot(morpho_plot, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=IUCN,
                 size=contributionPer), alpha=0.7) +
  scale_radius(name = "Contribution to FD (%)", range = c(0, 10), 
               breaks = c(0.005, 0.01, 0.02, 0.04, 0.06)) + theme_bw() + 
  scale_color_manual(name = "Period", 
                     values = rev(c("#D33F6A", "#4A6FE3", 
                                    "#BDC2E3"))) +
  labs(x='PC1 (73%)', y="PC2 (8%)") +
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14), 
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  ylim(c(-5, 4)) + xlim(c(-12, 7)) 


# Add density curves to y and x axis
xdens <- 
  axis_canvas(g1, axis = "x") + 
  geom_density(data = morpho_plot, trim = F, 
               aes(x = PC1, fill = IUCN), alpha = 0.7, color=NA) +
  scale_fill_manual(name = "", 
                    values = rev(c("#D33F6A", "#4A6FE3", 
                                   "#BDC2E3"))) +
  scale_color_manual(name = "", 
                     values = rev(c("#D33F6A", "#4A6FE3", 
                                    "#BDC2E3")))

ydens <-
  axis_canvas(g1, axis = "y", 
              coord_flip = TRUE) + 
  geom_density(data = morpho_plot, 
               aes(x = PC2, fill = IUCN), alpha = 0.7, color=NA) +
  scale_fill_manual(name = "", 
                    values = rev(c("#D33F6A", "#4A6FE3", 
                                   "#BDC2E3"))) +
  scale_color_manual(name = "", 
                     values = rev(c("#D33F6A", "#4A6FE3", 
                                    "#BDC2E3"))) +
  coord_flip() 

g2 <- g1 %>%
  insert_xaxis_grob(xdens, grid::unit(0.5, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(0.5, "in"), position = "right") %>%
  ggdraw()

jpeg("morphospace_density.jpeg", width = 18, height =12,
     res = 300, units = "cm")
g2
dev.off()

########################################
#####FIGURE 1 MAP#######################
########################################

g8b <- build_map(extinct = extinct)

jpeg(filename = "Map.jpeg", 
     width = 20, 
     height = 15, 
     units = "cm", 
     res = 300) # Save the plot out
g8b
dev.off()
