#### Meta data ####
# =============================================================================
# Project: Lake Fazilman Holocene climate, vegetation and human impact reconstructions
# Script: main.R
#
# Use:
# Change the Boolean variables (e.g., Import = T) to import, reconstruct past environmental
# variables and drow the figure related to the article (Dugerdil et al., 2026a)
# 
# Author: Lucas Dugerdil
# ORCID: 0000-0003-0266-564X
#
# Affiliations:
# 1) Univ. Lyon, ENS de Lyon, Université Lyon 1, CNRS, UMR 5276 LGL-TPE,
#    F-69364, Lyon, France
# 2) Université de Montpellier, CNRS, IRD, EPHE, UMR 5554 ISEM,
#    Montpellier, France
#
# Description:
# This script loads and cleans raw proxy data (XRF, magnetic susceptibility,
# spectrocolorimetry, pollen, NPPs, lipid biomarkers) from the Lake Fazilman
# sediment core to prepare datasets for subsequent analyses.
#
# Related publication:
# Dugerdil, L. (2026). Mid-Holocene wet optimum and Early-Late Holocene arid
# phases shaped steppe-forest vegetation and human societies of Uzbekistan:
# multi-proxy evidences from Lake Fazilman. CATENA. https://doi.org/10.1016/j.catena.2025.109759
#
# License:
# Creative Commons Attribution 4.0 International (CC BY 4.0)
#
# Citation:
# If you use this code, please cite:
# Lucas Dugerdil. (2025). LucasDugerdil/Paleo_Fazilman: v1.0.1 (v1.0.1). 
# Zenodo. https://doi.org/10.5281/zenodo.17381815
#
# Created: 2025-08-01
# Last modified: 2025-12-27
# =============================================================================

#### Import data ####
Import = F
if(Import == T){
  #### Scripts ####
  source("Import/Script/setup_packages.R")
  source("Import/Script/plot_functions.R")
  source("Import/Script/reconstruction_functions.R")
  
  #### Fazilman data ####
  ADM.Faz <- readRDS("Import/Data/Age_depth_model_Fazilman.Rds")
  XRF.Faz.MC <- readRDS("Import/Data/XRF_Fazilman.Rds")
  XRF.Faz.SD <- readRDS("Import/Data/XRF_Fazilman_SD.Rds")
  MS.Faz.MC <- readRDS("Import/Data/MS_Fazilman.Rds")
  Fazilman.PoF <- readRDS("Import/Data/Pollen_FA_Fazilman.Rds")
  Fazilman.NPPIx <- readRDS("Import/Data/NPP_influx_Fazilman.Rds")
  Fazilman.PoIx <- readRDS("Import/Data/Pollen_influx_Fazilman.Rds")
  Fazilman.AlF <- readRDS("Import/Data/Algue_FA_Fazilman.Rds")
  MP.Fazilman.Age <- readRDS("Import/Data/Pollen_age_Fazilman.Rds")
  Pol.Sum.Fazilman <- readRDS("Import/Data/Pollen_sum_Fazilman.Rds")
  GDGT.Faz <- readRDS("Import/Data/GDGT_Fazilman_Mean.Rds")
  GDGT.Faz.MaxI <- readRDS("Import/Data/GDGT_Fazilman_maxI.Rds")
  GDGT.Faz.MinI <- readRDS("Import/Data/GDGT_Fazilman_minI.Rds")
  Cores.metadata <- readRDS("Import/Data/Metadata_Fazilman.Rds")
  Units.L.GDGT <- readRDS("Import/Data/GDGT_Litho_Fazilman.Rds")
  
  Oversampling <- which(GDGT.Faz$Age < 2500)[which(GDGT.Faz$Age < 2500) %% 2 == 0]
  Keep.samples <- which(! seq(1, nrow(GDGT.Faz)) %in% Oversampling)
  GDGT.Faz <- GDGT.Faz[Keep.samples,]
  GDGT.Faz.MaxI <- GDGT.Faz.MaxI[Keep.samples,]
  GDGT.Faz.MinI <- GDGT.Faz.MinI[Keep.samples,]
  
  #### Other Paleo data ####
  HYDE.Uz <- readRDS("Import/Models/HYDE_uz_plot_1500.Rds")
  Datatime.MAP <- data.frame(read.csv("Import/Models/Uz_timeline_MAP.csv", sep = "\t", dec = ",", header = T))
  Datatime.MAAT <- data.frame(read.csv("Import/Models/Uz_timeline_MAAT.csv", sep = "\t", dec = ",", header = T))
  
  #### Modern calibration data ####
  Uz.TaxaCorresp  <- data.frame(read.csv(file="Import/Dictionary/CWM_table.csv",sep=",",dec=".", header=T, row.names=1))
  MT.Uz.full_gf <- readRDS("Import/Data_calibration/Modern_trait_data_Uzbekistan.Rds")
  PFT2biom <- readRDS("Import/Dictionary/PFT_2_Biome.Rds")
  Tax2PFT <- readRDS("Import/Dictionary/Taxon_2_PFT.Rds")
  
  TUSDB.MP <- readRDS("Import/Data_calibration/TUSDB_raw_data.Rds")
  ACADB.MP <- readRDS("Import/Data_calibration/ACADB_raw_data.Rds")
  Uz.MP.FT <- readRDS("Import/Data_calibration/TUSDB_raw_data.Rds")
  MAT.TUSDB <- readRDS("Import/Data_calibration/MAT_TUSDB.Rds")
  WAPLS.TUSDB <- readRDS("Import/Data_calibration/WAPLS_TUSDB.Rds")
  BRT.TUSDB <- readRDS("Import/Data_calibration/BRT_TUSDB.Rds")
  MAT.ACADB <- readRDS("Import/Data_calibration/MAT_ACADB.Rds")
  WAPLS.ACADB <- readRDS("Import/Data_calibration/WAPLS_ACADB.Rds")
  BRT.ACADB <- readRDS("Import/Data_calibration/BRT_ACADB.Rds")
  
  BRT.brACA <- readRDS("Import/Data_calibration/BRT_brACA_tuned.Rds")
  BRT.brACA.karid <- readRDS("Import/Data_calibration/BRT_brACA_karid_tuned.Rds")
  BRT.brACA.kwet <- readRDS("Import/Data_calibration/BRT_brACA_kwet_tuned.Rds")
  
  Msurf.Faz <- readRDS("Import/Data_calibration/Fazilman_200km_buff_BRT.Rds")
  Analog.Fazilman <- readRDS("Import/Data_calibration/Fazilman_200km_buff.Rds")
  Mclim.Faz <- readRDS("Import/Data_calibration/Fazilman_surf_clim.Rds")
  Meco.Faz <- readRDS("Import/Data_calibration/Fazilman_surf_eco.Rds")
  Cluster.prediction.ACADB.brGDGT <- readRDS("Import/Data_calibration/Cluster.prediction.ACADB.brGDGT.Rds")
  
  }

#### Reconstruction calculation ####
Reconstructions = F
if(Reconstructions == T){
  #### Pollen (CWM Trait reconstructions) ####
  CWM = F
  if(CWM == T){
    #### Clean pollen data for CWM calculation ####
    Fazilman.PoF_ss <- cbind(PT.ss = Uz.TaxaCorresp$PT_ss.label[match(row.names(Fazilman.PoF), row.names(Uz.TaxaCorresp))], Fazilman.PoF)
    Fazilman.PoF_ss <- aggregate(Fazilman.PoF_ss[-1], by = list(Fazilman.PoF_ss[["PT.ss"]]), sum)
    row.names(Fazilman.PoF_ss) <- Fazilman.PoF_ss$Group.1
    Fazilman.PoF_ss <- Fazilman.PoF_ss[-1]
    Fazilman.PoF_ss <- round(Fazilman.PoF_ss*100, digits = 2)
  
    #### CWM-traits calculation ####
    Seuil.pcover <- 70
    MT.Uz.MP_ss_gf <- MT.Uz.full_gf[MT.Uz.full_gf$Rank == "PT_ss", c(2,4:ncol(MT.Uz.full_gf))]
    MCWT.Faz.MP_ss_gf <- CWT.calculation(MT = MT.Uz.MP_ss_gf, MP = Fazilman.PoF_ss, Mclim = data.frame(t(MP.Fazilman.Age)), Accep.seuil = Seuil.pcover)
    saveRDS(MCWT.Faz.MP_ss_gf, "Results/CWM_Fazilman.Rds")
    }
  
  #### Pollen (Biomization) ####
  Biomization = F
  if(Biomization == T){
    PFT2biom <- PFT2biom[!PFT2biom$RECN %in% c("PION", "AQUA", "ANTH"),]
    row.names(Fazilman.PoF)[]
    MP_Fazilman.conv <- Fossil.MAT.prep(data.frame(t(Fazilman.PoF)), data.frame(t(MP.Fazilman.Age)), Type_MAT = "Type_Odile",  Corresp_name = Uz.TaxaCorresp, Displot = F, Show.message = T)
    
    Fazilman.biomiz <- Biomization_cal(Pollen = MP_Fazilman.conv, Tax2PFT = Tax2PFT, PFT2biom = PFT2biom,
                                 Display.warning = T, Verbose = T, Posterior.PFT.correction = F, bioindic = F,
                                 Save.path = "Results/Biomes_Fazilman.Rds")
    Fazilman.biomiz <- Fazilman.biomiz$Biom.output
    
  }
  
  #### Pollen (Climate reconstructions) ####
  Climate = F
  if(Climate == T){
    Fazilman.MA <- data.frame(t(MP.Fazilman.Age))
    Fazilman.MA <- Fazilman.MA[row.names(Fazilman.MA) %in% names(Fazilman.PoF),]
    MP_Fazilman.conv <- Fossil.MAT.prep(data.frame(t(Fazilman.PoF)), Fazilman.MA, Type_MAT = "Type_Odile",  Corresp_name = Uz.TaxaCorresp, Show.message = F)
    
    #### TUSDB calibration (Uzbekistan + Tajikistan modern pollen assemblages) ####
    TUSDB = F
    if(TUSDB == T){
      MP_Fazilman.conv.TUSDB <- Fossil.corresp.surface(MP_Fazilman.conv, TUSDB.MP)
      Fazilman.TUSDB <- FT.core(
        Model.WAPLS = WAPLS.TUSDB,
        Model.MAT = MAT.TUSDB,
        Model.BRT = BRT.TUSDB,
        MCore = MP_Fazilman.conv.TUSDB,
        MAge = Fazilman.MA$Age,
        Fit.val = 0.25, LakeName = "Fazilman",
        Only.fit = F, Model.param.show = T, Save.RDS = T, Displot = F,
        Save.path = "Results/Fazilman.csv")}
    
    #### ACADB calibration (Arid Central Asian modern pollen assemblages) ####
    ACADB = F
    if(ACADB == T){
      MP_Fazilman.conv.ACADB <- Fossil.corresp.surface(MP_Fazilman.conv, ACADB.MP)
      Fazilman.ACADB <- FT.core(
        Model.WAPLS = WAPLS.ACADB,
        Model.MAT = MAT.ACADB,
        Model.BRT = BRT.ACADB,
        MCore = MP_Fazilman.conv.ACADB,
        MAge = Fazilman.MA$Age,
        Fit.val = 0.25, LakeName = "Fazilman",
        Only.fit = F, Model.param.show = T, Save.RDS = T, Displot = F,
        Save.path = "Results/Fazilman.csv")}
    
  }
  
  #### brGDGT (Boosted Regression Trees calibrations + Local scaling) ####
  brGDGT = F
  if(brGDGT == T){
    #### Clean data GDGTs ####
    names(GDGT.Faz)[names(GDGT.Faz) == "Top"] <- "MCD"
    names(GDGT.Faz.MaxI)[names(GDGT.Faz.MaxI) == "Top"] <- "MCD"
    names(GDGT.Faz.MinI)[names(GDGT.Faz.MinI) == "Top"] <- "MCD"
    GDGT.Faz$MCD <- 10*GDGT.Faz$MCD
    GDGT.Faz.MaxI$MCD <- 10*GDGT.Faz.MaxI$MCD
    GDGT.Faz.MinI$MCD <- 10*GDGT.Faz.MinI$MCD
    
    GDGT.Faz$pCren <- 100*GDGT.Faz$pCren
    GDGT.Faz.MaxI$pCren <- 100*GDGT.Faz.MaxI$pCren
    GDGT.Faz.MinI$pCren <- 100*GDGT.Faz.MinI$pCren
    
    GDGT.Faz$Litho <- Units.L.GDGT
    GDGT.Faz.MaxI$Litho <- Units.L.GDGT
    GDGT.Faz.MinI$Litho <- Units.L.GDGT
    
    #### Add ternary information ####
    GDGT.Faz.tern <- Cal.Mtern(GDGT.Faz)
    GDGT.Faz <- cbind(GDGT.Faz, GDGT.Faz.tern)
    GDGT.Faz.tern <- Cal.Mtern(GDGT.Faz.MaxI)
    GDGT.Faz.MaxI <- cbind(GDGT.Faz.MaxI, GDGT.Faz.tern)
    GDGT.Faz.tern <- Cal.Mtern(GDGT.Faz.MinI)
    GDGT.Faz.MinI <- cbind(GDGT.Faz.MinI, GDGT.Faz.tern)
    
    #### Machine-Learning reconstructions (Boosted Regression Trees) ####
    GDGT.Faz.conv <- GDGT.Faz[grep("^f.I", names(GDGT.Faz))]
    GDGT.Faz.conv <- GDGT.Faz.conv[!grepl("7Me", names(GDGT.Faz.conv))]
    
    Fazilman.brACA <- FT.core(Model.BRT = BRT.brACA,
                              MCore = GDGT.Faz.conv, MAge = GDGT.Faz$Age, GDGT.model = "BRT",
                              LakeName = "Fazilman", Only.fit = T, Save.RDS = T, Displot = F, GDGT = T,
                              Save.path = "Results/Fazilman.csv")
    
    Fazilman.brACA.karid <- FT.core(Model.BRT = BRT.brACA.karid,
                                    MCore = GDGT.Faz.conv, MAge = GDGT.Faz$Age, GDGT.model = "BRT",
                                    LakeName = "Fazilman", Only.fit = T, Save.RDS = T, Displot = F, GDGT = T,
                                    Save.path = "Results/Fazilman.csv")
    
    Fazilman.brACA.kwet <- FT.core(Model.BRT = BRT.brACA.kwet,
                                   MCore = GDGT.Faz.conv, MAge = GDGT.Faz$Age, GDGT.model = "BRT",
                                   LakeName = "Fazilman", Only.fit = T, Save.RDS = T, Displot = F, GDGT = T,
                                   Save.path = "Results/Fazilman.csv")
    Fazilman.BRT.ACADB <- readRDS("Results/Fazilman_BRT_brACA.Rds")[[2]]
    Fazilman.BRT.karid <- readRDS("Results/Fazilman_BRT_karid.Rds")[[2]]
    Fazilman.BRT.kwet  <- readRDS("Results/Fazilman_BRT_kwet.Rds")[[2]]
    
    #### Add FTs results on brGDGTs ####
    GDGT.Faz.ML <- Add.ML.to.cores(
      Mcore = GDGT.Faz,
      Mmin = GDGT.Faz.MinI,
      Mmax = GDGT.Faz.MaxI,
      Mcore.ML = list(
        "BRT_ACA" = Fazilman.BRT.ACADB, 
        "BRT_KAR" = Fazilman.BRT.karid,
        "BRT_KWT" = Fazilman.BRT.kwet
      ),
      
      Msurf.ML = list(
        "BRT_ACA" = BRT.brACA, 
        "BRT_KAR" = BRT.brACA.karid,
        "BRT_KWT" = BRT.brACA.kwet
      ),
    )
    
    GDGT.Faz <- GDGT.Faz.ML[[1]]
    GDGT.Faz.MinI <- GDGT.Faz.ML[[2]]
    GDGT.Faz.MaxI <- GDGT.Faz.ML[[3]]
    
    #### Local rescale ####
    Add.LR = T
    if(Add.LR == T){
      #### Local correction ####
      Rescale.Faz <- GDGT.local.rescale(List.params = c("MAAT", "AI", "MPCOQ", "MAF"),
                                        Keep.models = c("BRT"),
                                        Msurf = Msurf.Faz,
                                        Mclim = Mclim.Faz,
                                        Mpaleo = GDGT.Faz,
                                        Mpaleo.MaxI = GDGT.Faz.MaxI,
                                        Mpaleo.MinI = GDGT.Faz.MinI)
      Msurf.Faz <- Rescale.Faz[[1]]
      GDGT.Faz <- Rescale.Faz[[2]]
      GDGT.Faz.MaxI <- Rescale.Faz[[3]]
      GDGT.Faz.MinI <- Rescale.Faz[[4]]
      
      #### Cluster prediction (surf) ####
      Cm.predict.Faz = T
      if(Cm.predict.Faz == T){
        M1 <- setNames(Msurf.Faz[grepl("BRT_ACA", names(Msurf.Faz))&!grepl("_LS", names(Msurf.Faz))], gsub("_BRT_ACA", "", names(Msurf.Faz[grepl("BRT_ACA", names(Msurf.Faz))&!grepl("_LS", names(Msurf.Faz))])))
        M2 <- setNames(Msurf.Faz[grepl("BRT_KAR", names(Msurf.Faz))&!grepl("_LS", names(Msurf.Faz))], gsub("_BRT_KAR", "", names(Msurf.Faz[grepl("BRT_KAR", names(Msurf.Faz))&!grepl("_LS", names(Msurf.Faz))])))
        M3 <- setNames(Msurf.Faz[grepl("BRT_KWT", names(Msurf.Faz))&!grepl("_LS", names(Msurf.Faz))], gsub("_BRT_KWT", "", names(Msurf.Faz[grepl("BRT_KWT", names(Msurf.Faz))&!grepl("_LS", names(Msurf.Faz))])))
        Save.BRT <- "Results/Faz_BRT_CB_surf_buffer.Rds"
        Met1 <- "BRT"
        
        M4 <- setNames(Msurf.Faz[grepl("BRT_ACA_LS", names(Msurf.Faz))], gsub("_BRT_ACA_LS", "", names(Msurf.Faz[grepl("BRT_ACA_LS", names(Msurf.Faz))|grepl("Age", names(Msurf.Faz))])))
        M5 <- setNames(Msurf.Faz[grepl("BRT_KAR_LS", names(Msurf.Faz))], gsub("_BRT_KAR_LS", "", names(Msurf.Faz[grepl("BRT_KAR_LS", names(Msurf.Faz))|grepl("Age", names(Msurf.Faz))])))
        M6 <- setNames(Msurf.Faz[grepl("BRT_KWT_LS", names(Msurf.Faz))], gsub("_BRT_KWT_LS", "", names(Msurf.Faz[grepl("BRT_KWT_LS", names(Msurf.Faz))|grepl("Age", names(Msurf.Faz))])))
        Save.BRT.LS <- "Results/Faz_BRT_CB_LS_surf_buffer.Rds"
        Met2 <- "BRT-LS"
        
        Param.clim.order <- c("MPCOQ", "MAAT", "MAF", "AI")
        
        M1 <-M1[Param.clim.order]
        M2 <-M2[Param.clim.order]
        M3 <-M3[Param.clim.order]
        M4 <-M4[Param.clim.order]
        M5 <-M5[Param.clim.order]
        M6 <-M6[Param.clim.order]
        
        M1 <- cbind(ID = seq(1: nrow(M1)), M1)
        M2 <- cbind(ID = seq(1: nrow(M2)), M2)
        M3 <- cbind(ID = seq(1: nrow(M3)), M3)
        M4 <- cbind(ID = seq(1: nrow(M4)), M4)
        M5 <- cbind(ID = seq(1: nrow(M5)), M5)
        M6 <- cbind(ID = seq(1: nrow(M6)), M6)
        Msurf.Faz <- cbind(ID = seq(1: nrow(Msurf.Faz)), Msurf.Faz)
        
        #### Plot BRT ####
        Faz.surf.BRT.CB <- Combine.ML.cluster(
          List.models = list(M1 = M1, M2 = M2, M3 = M3),
          Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
          Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
          GDGT.paleo = Msurf.Faz, Method = Met1,
          Surf.val = Cores.metadata$MAAT[row.names(Cores.metadata) == "Fazilman"],
          Compare.curve = c("MAAT_mr_DJ", "MAAT_DJ_5Me", "MAAT_Chen_Tjk"), #Core.name = "Fazilman",
          Plot.y = "ID", Plot.y.lab = "Age (yr cal BP)", Param.clim = "MAAT", Cluster.prob = "K-cold/wet",
          Facet = F, Only.best = F, Show.proba = F, Highlight.combined = F,
          Save.path = Save.BRT)
        
        #### Plot BRT-LS ####
        Faz.surf.BRT.CB.LS <- Combine.ML.cluster(
          List.models = list(M1 = M4, M2 = M5, M3 = M6),
          Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
          Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
          GDGT.paleo = Msurf.Faz, Method = Met2,
          Surf.val = Cores.metadata$MAAT[row.names(Cores.metadata) == "Fazilman"],
          Compare.curve = c("MAAT_mr_DJ", "MAAT_DJ_5Me", "MAAT_Chen_Tjk"), Core.name = "Fazilman",
          Plot.y = "ID", Plot.y.lab = "Age (yr cal BP)", Param.clim = "MAAT", Cluster.prob = "K-cold/wet",
          Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T, 
          Save.path = Save.BRT.LS)
      }
      
      Faz.surf.BRT.CB <- readRDS("Results/Faz_BRT_CB_surf_buffer.Rds")
      Faz.surf.BRT.CB.LS <- readRDS("Results/Faz_BRT_CB_LS_surf_buffer.Rds")
      
      #### Add CB et CB LS to surf ####
      Msurf.Faz <- Add.ML.to.cores(
        Mcore = Msurf.Faz,
        Mcore.ML = list(
          "BRT_ACA_CB" = Faz.surf.BRT.CB[[1]],
          "BRT_ACA_CB_LS" = Faz.surf.BRT.CB.LS[[1]]
        ),
      )
      Msurf.Faz <- Msurf.Faz[[1]]
    }
    else{
      # Msurf.Faz <- readRDS("Results/Fazilman_200km_buff_BRT_LS_CB.Rds")
      # Mclim.Faz <- readRDS("Results/Fazilman_surf_clim.Rds")
      # Meco.Faz  <- readRDS("Results/Fazilman_surf_eco.Rds")
      
      # GDGT.Faz <- readRDS("Results/Fazilman_GDGT_full_models_LS.Rds")
      # GDGT.Faz.MaxI <- readRDS("Results/Fazilman_GDGT_full_models_MaxI_LS.Rds")
      # GDGT.Faz.MinI <- readRDS("Results/Fazilman_GDGT_full_models_MinI_LS.Rds")
    }
    
    #### Cluster prediction ####
    Cm.predict.Faz = T
    if(Cm.predict.Faz == T){
      Use.LS = T
      #### Use LS or not ####
      M1 <- Fazilman.BRT.ACADB
      M2 <- Fazilman.BRT.karid
      M3 <- Fazilman.BRT.kwet
      Plot.BRT <- "Figures/Figure_S9.pdf"
      Save.BRT <- "Results/Faz_BRT_CB.Rds"
      Met1 <- "BRT"
      M4 <- setNames(GDGT.Faz[grepl("BRT_ACA_LS", names(GDGT.Faz))|grepl("Age", names(GDGT.Faz))], gsub("_BRT_ACA_LS", "", names(GDGT.Faz[grepl("BRT_ACA_LS", names(GDGT.Faz))|grepl("Age", names(GDGT.Faz))])))
      M5 <- setNames(GDGT.Faz[grepl("BRT_KAR_LS", names(GDGT.Faz))|grepl("Age", names(GDGT.Faz))], gsub("_BRT_KAR_LS", "", names(GDGT.Faz[grepl("BRT_KAR_LS", names(GDGT.Faz))|grepl("Age", names(GDGT.Faz))])))
      M6 <- setNames(GDGT.Faz[grepl("BRT_KWT_LS", names(GDGT.Faz))|grepl("Age", names(GDGT.Faz))], gsub("_BRT_KWT_LS", "", names(GDGT.Faz[grepl("BRT_KWT_LS", names(GDGT.Faz))|grepl("Age", names(GDGT.Faz))])))
      # Plot.BRT.LS <- "Figures/BRT_comb_Faz_LS.pdf"; 
      Plot.param <- "MAAT"
      Save.BRT.LS <- "Results/Faz_BRT_CB_LS.Rds"
      Met2 <- "BRT-LS"
      
      Param.clim.order <- c("Age", "MPCOQ", "MAAT", "MAF", "AI")
      M1 <-M1[Param.clim.order]
      M2 <-M2[Param.clim.order]
      M3 <-M3[Param.clim.order]
      M4 <-M4[Param.clim.order]
      M5 <-M5[Param.clim.order]
      M6 <-M6[Param.clim.order]
      
      #### Plot BRT ####
      Faz.BRT.CB <- Combine.ML.cluster(
        List.models = list(M1 = M1, M2 = M2, M3 = M3),
        Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
        Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
        GDGT.paleo = GDGT.Faz, 
        Surf.val = Cores.metadata$MAAT[row.names(Cores.metadata) == "Fazilman"],
        Compare.curve = c("MAAT_mr_DJ", "MAAT_DJ_5Me", "MAAT_Chen_Tjk"),
        Core.name = "Fazilman",
        Plot.y = "Age", Plot.y.lab = "Time (cal. year BP)", Param.clim = Plot.param,
        Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T,
        Save.path = Save.BRT,
        H = 1100, W = 450, Save.plot = Plot.BRT)
      
      #### Plot BRT-LS ####
      Faz.BRT.CB.LS <- Combine.ML.cluster(
        List.models = list(M1 = M4, M2 = M5, M3 = M6),
        Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
        Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
        GDGT.paleo = GDGT.Faz, #Time.lim = c(0, 2700),
        Surf.val = Cores.metadata$MAAT[row.names(Cores.metadata) == "Fazilman"],
        Compare.curve = c("MAAT_mr_DJ", "MAAT_DJ_5Me", "MAAT_Chen_Tjk"), Core.name = "Fazilman",
        Plot.y = "Age", Plot.y.lab = "Age (yr cal BP)", Param.clim = "MAAT", Cluster.prob = "K-warm/arid",
        Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T,
        Save.path = Save.BRT.LS)
      
    }
    Faz.BRT.CB <- readRDS("Results/Faz_BRT_CB.Rds")
    Faz.BRT.CB.LS <- readRDS("Results/Faz_BRT_CB_LS.Rds")
    
    #### Determine RMSE ####
    RMSE.BRT.CB    <- RMSE.comb(ML.comb = Faz.BRT.CB, BRT.K1 = BRT.brACA.karid, BRT.K2 = BRT.brACA.kwet)
    RMSE.BRT.CB.LS <- RMSE.comb(ML.comb = Faz.BRT.CB.LS, BRT.K1 = BRT.brACA.karid, BRT.K2 = BRT.brACA.kwet)
    
    #### Add FTs combined on brGDGTs ####
    GDGT.Faz.ML <- Add.ML.to.cores(
      Mcore = GDGT.Faz,
      Mmin = GDGT.Faz.MinI,
      Mmax = GDGT.Faz.MaxI,
      Mcore.ML = list(
        "BRT_ACA_CB_LS" = Faz.BRT.CB.LS[[1]], 
        "BRT_ACA_CB" = Faz.BRT.CB[[1]]
      ),
      
      Msurf.ML = list(
        "BRT_ACA_CB_LS" = RMSE.BRT.CB.LS,
        "BRT_ACA_CB" = RMSE.BRT.CB
      ),
    )
    
    GDGT.Faz <- GDGT.Faz.ML[[1]]
    GDGT.Faz.MinI <- GDGT.Faz.ML[[2]]
    GDGT.Faz.MaxI <- GDGT.Faz.ML[[3]]
    
    #### Ensemble model test ####
    Ensemble.models = T
    if(Ensemble.models == T){
      ENS.Faz <- GDGT.ensemble.RMSE(
        Mpaleo = GDGT.Faz, Mpaleo.min = GDGT.Faz.MinI, Mpaleo.max = GDGT.Faz.MaxI,
        Msurf = Msurf.Faz, Mclim = Mclim.Faz, Ensemble.BAY = F, Ensemble.WM = T, Variance.method = "quadratic",
        Param.clim = c("MAAT", "AI",  "MAF", "MPCOQ"),
        Ensemble = list(c("MAAT_Chen_Tjk", "MAAT_mr_DJ", "MAAT_BRT_ACA_LS", "MAAT_BRT_ACA_CB_LS"),
                        c("AI_BRT_ACA_LS", "AI_BRT_ACA_CB_LS", "AI_BRT_KAR_LS"),
                        c("MAF_MSosa", "MAF_BRT_ACA_LS", "MAF_BRT_ACA_CB_LS", "MAF_full_Raberg"),
                        c("MPCOQ_BRT_ACA_LS", "MPCOQ_BRT_ACA_CB_LS")),
        Save.path = "Results/Faz_GDGT_full_models_CB_LS_ens.Rds"
      )}
    
    #### randomTF (GDGT) ####
    test.randomTF = F
    if(test.randomTF == T){
      Faz.rTF.ACADB.brGDGT <- Plot.randomTF(MPsurf = M.br.GDGT, MPpaleo = GDGT.Faz.conv,  Mclim = MC, Database = "ACADB", Lake = "Fazilman", Plot.MAT = F, Plot.WAPLS = F, 
                                            Save.path = "Resultats/Uzbekistan/GDGT/Cores/Fazilman/Func_trans/Fazilman_randomFT_ACADB.Rds", return.plot = T, Plot.RF = T, Plot.BRT = T,
                                            H = 500, W = 2000, Save.plot = "Figures/Uzbekistan/GDGT/Fazilman/randomTF/Fazilman_randomTF_ACADB.pdf")
      saveRDS(Faz.rTF.ACADB.brGDGT, "Resultats/Uzbekistan/GDGT/Cores/Fazilman/Func_trans/Fazilman_randomFT_ACADB_brGDGT_plot.Rds")
      
      Faz.rTF.ACADB.comb <- Plot.randomTF(MPsurf = M.br.GDGT.kwet, MPpaleo = GDGT.Faz.conv,  Mclim = MC.kwet, Database = "ACADB-combined", Lake = "Fazilman", Plot.MAT = F, Plot.WAPLS = F, 
                                          Save.path = "Resultats/Uzbekistan/GDGT/Cores/Fazilman/Func_trans/Fazilman_randomFT_ACADB_comb.Rds", return.plot = T, Plot.RF = T, Plot.BRT = T,
                                          H = 500, W = 2000, Save.plot = "Figures/Uzbekistan/GDGT/Fazilman/randomTF/Fazilman_randomTF_ACADB_comb.pdf")
      saveRDS(Faz.rTF.ACADB.comb, "Resultats/Uzbekistan/GDGT/Cores/Fazilman/Func_trans/Fazilman_randomFT_ACADB_comb_plot.Rds")
    }
    
    #### randomTF (merge with pollen) #### 
    Full.table = F
    if(Full.table == T){
      Faz.rTF.ACADB <- readRDS("Resultats/Uzbekistan/GDGT/Cores/Fazilman/Func_trans/Fazilman_randomFT_ACADB.Rds")
      Faz.rTF.ACADB.comb <- readRDS("Resultats/Uzbekistan/GDGT/Cores/Fazilman/Func_trans/Fazilman_randomFT_ACADB_comb.Rds")
      
      Faz.rTF.full <- rbind(Faz.rTF.ACADB, Faz.rTF.ACADB.comb) 
      
      Position.clim <- match(unique(Faz.rTF.full$Database), Faz.rTF.full$Database)
      Faz.rTF.full$Database <- as.character(Faz.rTF.full$Database)
      # Faz.rTF.full$Database[setdiff(seq(1:nrow(Faz.rTF.full)),Position.clim)] <- ""
      names(Faz.rTF.full)[ncol(Faz.rTF.full)] <- "$\\mathrm{PCA_{1} (var. \\%)^{(b)}}$"
      names(Faz.rTF.full)[ncol(Faz.rTF.full)-1] <- "Threshold$^{(a)}$"
      Faz.rTF.full <- Faz.rTF.full[Faz.rTF.full$Model == "BRT",]
      saveRDS(Faz.rTF.full, "Resultats/Uzbekistan/GDGT/Cores/Fazilman/Func_trans/Fazilman_randomFT_all_brGDGT.Rds")
    }
    
  }
  
}
if(Reconstructions == F){
  MCWT.Faz.MP_ss_gf <- readRDS("Results/CWM_Fazilman.Rds")
  Fazilman.biomiz <- readRDS("Results/Biomes_Fazilman_biomization.Rds")
  
  GDGT.Fazilman <- readRDS("Results/Faz_GDGT_full_models_CB_LS_ens.Rds")
  GDGT.Fazilman.MaxI <- readRDS("Results/Faz_GDGT_full_models_CB_LS_ens_MaxI.Rds")
  GDGT.Fazilman.MinI <- readRDS("Results/Faz_GDGT_full_models_CB_LS_ens_MinI.Rds")
  
  if(exists("Fazilman.ACADB") == F){Fazilman.ACADB <- readRDS("Results/Fazilman_ACADB.Rds")}
  Fazilman.TUSDB <- readRDS("Results/Fazilman_TUSDB.Rds")
  }

#### Figures plotting ####
Figures = T
if(Figures == T){
  #### Figure 3 ####
  Fig3 = F
  if(Fig3 == T){
    XRF.Faz.MC.onlyPCA <- subset(XRF.Faz.MC, select = c(MCD, Age, PC1, PC2))
    MS.Faz.MC.onlyMS <- subset(MS.Faz.MC, select = c(MCD, Age, MS))
    CONISS.Faz <- Plot.geoch.core(XRF = XRF.Faz.MC.onlyPCA, MS = MS.Faz.MC.onlyMS,
                                  Plot.x ="MCD",
                                  CONISS = T, Nzone = 8,
                                  Print.cluster.zone = T, Groupes = list(c("PC1", "MS"), c("PC2")),
                                  Select.interv.x = 50)

    CONISS.Faz$Min.Age <- ADM.Faz$mean[ADM.Faz$MCD %in% CONISS.Faz$Min.CONISS]
    CONISS.Faz$Max.Age <- ADM.Faz$mean[ADM.Faz$MCD %in% CONISS.Faz$Max.CONISS]
    
    Figure.3.plot <- Plot.geoch.core(XRF = XRF.Faz.MC, XRF.SD = XRF.Faz.SD,
                                     MS = MS.Faz.MC,
                                     GDGT = GDGT.Fazilman, GDGT.min = GDGT.Fazilman.MinI, GDGT.max = GDGT.Fazilman.MaxI,
                                     Manual.line = c(CONISS.Faz$Min.CONISS, max(CONISS.Faz$Max.CONISS)),
                                     Select.interv.x = 200, Plot.x ="MCD", Nb.ticks.minors = 25,
                                     Groupes = list(
                                       OM = c("Br/Ti", "Inc/Coh", "S/Ti", "RABD660_670"),
                                       Authigenic.Lithogenic = c("PC1", "MS"),
                                       Carbonates = c("Ca/Ti", "L"),
                                       Sand = c("Si/Ti", "Zr/Ti"),
                                       Bio.carbonates = c("Ca/Mg"),
                                       Salinity = c("IR7Me", "Sr/Ca"),
                                       Oxidation = c("Fe/Mn", "Q7.4"),
                                       Water.depth = c("pCren"),
                                       Anaeroby = c("S/Fe"),
                                       Alkalinity = c("IR6Me"),
                                       Soil.influx = c("PC2", "Rb/Sr", "BIT", "Ib.Ia")
                                     ),
                                     Manual.color.curves = c(
                                       Authigenic.Lithogenic = "darkred",
                                       Soil.influx = "#825c28",
                                       Salinity = "#C2294A",
                                       OM = "#468237",
                                       Anaeroby = "#9E0142",
                                       Oxidation = "#DF4D4B",
                                       Sand = "#C3B783",
                                       Bio.carbonates = "#93D3A4",
                                       Carbonates = "#5E4FA2",
                                       Water.depth = "#4075B4",
                                       Alkalinity = "#439BB5"),
                                     Keep.MS = c("L", "MS", "RABD660_670", "Q7.4"),
                                     Ratio = c("Sr/Ca", "Br/Ti", "Inc/Coh", "S/Fe","Fe/Mn", "Si/Ti", "Zr/Ti", "Ca/Mg", "Ca/Ti"),
                                     Keep.XRF = c("PC1", "PC2"),
                                     Keep.GDGT = c("BIT", "Ib.Ia", "pCren", "IR7Me", "IR6Me"),
                                     Arrows = T,  Reverse.arrow = c(3,14), Arrow.y = 0.032,
                                     Temp.zone = rev(Fazilman.clus.color), Col.area = 1,
                                     Zone.clim = sort(c(CONISS.Faz$Min.CONISS, CONISS.Faz$Max.CONISS)),
                                     Name.zone = CONISS.Faz$hclust_zone,
                                     C14.path = "Import/Data/C14_Fazilman.csv",
                                     H = 700/1.2, W = 1900/1.2, Save.plot = "Figures/Figure_3.pdf")
    
  }
  
  #### Figure 4 ####
  Fig4 = F
  if(Fig4 == T){
    #### Build A ####
    Figure.4A.plot <- DP_pol(Fazilman.PoF, MP.Fazilman.Age, Plot.x = "Age",
                             Mpol_conc = Fazilman.PoIx,
                             Malg_frac = Fazilman.AlF,
                             NPP_lake = Fazilman.NPPIx,
                             Path.index = "Import/Dictionary/Pollen_table.csv",
                             Auto.Pollen.lab = T, Al.Pol.influx = F,
                             Zone.clim = c(50, 550, 650, 950, 1600, 1800, 1900, 2500, 3750, 4400, 5000, 7500, 8200, max(MP.Fazilman.Age[3,])),
                             Temp.zone = c("C","W","C","W", "C", "W", "C"),
                             Name.zone = c("LIA", "WMP", "DACP", "RWP", "4.2k", "Mid-Hol.", "Early Hol."),
                             Seuil.Pour = 0.03, Seuil.Pour.alg = 0.06,
                             Time.window = 500,
                             Seuil.NPP = 275, CONISS = T, Nzone = 6, Show.stats = F,
                             Order.algues =  c(8,5,4,3,6,1,11,10,2,7),
                             Order.pollen = c(5,1,2,3,4,6,7:15), 
                             Order.fungal = c(2:9,1,12,11,10),
                             Keep.cortege = c("Grazzing_Spore", "Saprophilous_Spore", "Aridity_Spore"),
                             AP.NAP = F, Only.AP = F, Sort.taxon = "Auto", Pollen.diversity = F, Spore.diversity = T,
                             Pol.influx.total = T, Spore.influx = F, Cortege = T, Richesse.type = T, Ratio = TabRa,
                             Save.plot = "Figures/Figure_4A.pdf", H = 1100, W = 3200)
    
    #### Build B ####
    B <- DP_plot_Fungal(Fazilman.NPPIx, MP.Fazilman.Age, Conc.Plot = F,
                        Pol.Sum = Pol.Sum.Fazilman, Auto.lab.NPP = T, Time.window = 500,  
                        Cortege = T, Col.cortege = T, Log.trans = F,
                        Limites = c(-100, 1400),
                        Path.index = "Import/Dictionary/NPP_table.csv",
                        CONISS = F, nlake = "Fazilman", Nzone = 6, Seuil.NPP = 1, Sort.taxon = "Auto.mean",
                        Spore.influx = T, Spore.diversity = T, Diversity.pollen.ratio = F, Richesse.type = T, 
    )
    
    #### Extraction and Clean de B #### 
    Mgg <- cbind(B[[1]], Age = t(B[[2]]))
    Mgg <- reshape2::melt(Mgg, id = c("Age"))
    Mgg$variable <- gsub("cf\\.", "cf", Mgg$variable)
    Mgg <- Mgg[Mgg$value != 0,]
    Mgg$Labels <- B[[4]]$Label[match(Mgg$variable, B[[4]]$Nom)]
    Mgg$Couleur <- B[[4]]$Couleur[match(Mgg$variable, B[[4]]$Nom)]
    Mgg$Cortege <- B[[4]]$Cortege[match(Mgg$variable, B[[4]]$Nom)]
    Mgg$variable[Mgg$variable == "Parasitic_Spore"] <- "Plant pathogens"
    Mgg$Labels[Mgg$Labels == "Parasitic spores"] <- "Plant pathogens"
    
    Mgg$Cortege[Mgg$Cortege == "NaN"] <- "Other NPPs"
    Mgg <- Mgg[order(Mgg$Cortege, Mgg$Labels),]
    Mgg$variable <- factor(Mgg$variable, unique(Mgg$variable), ordered = T)
    Or2 <- unique(Mgg[c("Labels", "Cortege")])
    Or2$Bool <- F 
    Or2$Bool[Or2$Labels == Or2$Cortege] <- T 
    Or2 <- Or2[order(Or2$Bool),]  
    Or2 <- Or2[c(1:(grep("influx", Or2$Labels)-1),(grep("influx", Or2$Labels)+1):nrow(Or2), grep("influx", Or2$Labels)),]
    Mgg$Labels <- factor(Mgg$Labels, Or2$Labels, ordered = T)
    Mgg$value <- as.numeric(Mgg$value)
    Mgg$Age <- as.numeric(Mgg$Age)
    My_lab <- unique(Mgg$Cortege)
    My_color <- c(unique(Mgg$Couleur), "grey")
    names(My_color) <- unique(Mgg$Cortege)
    Time.window = 100
    Limites = c(-100, 1400)
    X.lab.size = 8
    Display.legend = "none"
    Line.coniss <- geom_hline(aes(yintercept = 160))
    
    #### Plot histo fungal ####
    Plot.compa <- ggplot(Mgg, aes(x = value, y = Age, fill = Cortege)) +
      # Line.coniss +
      geom_vline(xintercept = 0, color = "grey70", size = .2)+
      geom_colh(width = 25, position = "dodgev", na.rm = T) +
      facet_abundanceh(vars(Labels), rotate_facet_labels = 45,
                       dont_italicize = c("spp?\\.", "type", "[Oo]ther", "^(.*?)eae", ".* spores", ".* markers", ".* pathogens", "^HdV.*", "cf.", "undet.")) +
      scale_fill_manual(values = My_color, labels = My_lab)+
      scale_y_reverse(breaks = seq(Limites[1],Limites[2], Time.window))+ ylab("Age (cal. year BP)")+
      xlab(expression(10^3~"x"~"#"*grains*.*cm^-3))+
      scale_x_continuous(breaks = round(seq(0, max(Mgg$value), length.out = 5), digits = -2), expand = c(0,0))+
      theme(panel.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), strip.background = element_blank(),
            axis.text.x = element_text(size = X.lab.size, angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 13), strip.text = element_text(size = 9),
            axis.line.x = element_line(lineend = "butt", color = "grey70", linewidth = .1),
            axis.ticks.x = element_line(lineend = "butt", color = "grey70", linewidth = .1),
            legend.position = Display.legend, legend.background = element_blank(), legend.key = element_blank(), legend.box.background = element_blank(),
            panel.spacing = unit(.9, "lines"))
    
    #### Clean data diversity ####
    Mgg2 <- cbind(B[[3]], Age = t(B[[2]]))
    Mgg2 <- Mgg2[-c(1)]
    
    Mgg2 <- reshape2::melt(Mgg2, id = c("Age"))
    Mgg2$Couleur <- B[[4]]$Couleur[match(gsub("\\.", " ", Mgg2$variable), B[[4]]$Nom)]
    Mgg2$Cortege <- B[[4]]$Cortege[match(Mgg2$variable, B[[4]]$Nom)]
    levels(Mgg2$variable)[levels(Mgg2$variable) == "Peat_Spore"] <- "Wetland markers"
    levels(Mgg2$variable)[levels(Mgg2$variable) == "Erosion_Spore"] <- "Erosion markers"
    levels(Mgg2$variable)[levels(Mgg2$variable) == "Aridity_Spore"] <- "Aridity markers"
    levels(Mgg2$variable)[levels(Mgg2$variable) == "Parasitic_Spore"] <- "Plant pathogens"
    levels(Mgg2$variable)[levels(Mgg2$variable) == "Saprophilous_Spore"] <- "Saprophytic spores"
    levels(Mgg2$variable)[levels(Mgg2$variable) == "Grazzing_Spore"] <- "Dung fungal spores"
    levels(Mgg2$variable)[levels(Mgg2$variable) == "Other.NPPs"] <- "Other NPPs"
    Or3 <- c("Aridity markers", "Dung fungal spores", "Erosion markers", "Plant pathogens", "Saprophytic spores", "Wetland markers", "Other NPPs")
    Mgg2$variable <- factor(Mgg2$variable, rev(Or3), ordered = T)
    
    #### Plot spore diversity ####
    Plot.div <- ggplot(Mgg2, aes(x = value, y = Age, fill = factor(variable), group = variable)) +
      scale_y_reverse(breaks = seq(Limites[1],Limites[2], Time.window))+
      geom_colh(width = 25, color = NA, fill = NA, na.rm = T) +
      scale_fill_manual(values = My_color, name = "Indicator groups")+
      geom_area(data = Mgg2, mapping = aes(x = value, y = Age, fill = factor(variable), group = factor(variable)), position = "stack", orientation = "y")+
      xlab("# Taxa")+
      theme(panel.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), strip.background = element_blank(),
            axis.text.x = element_text(size = X.lab.size, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_blank(), axis.title.y = element_blank(), #axis.ticks.y = element_blank(),
            axis.line = element_line(lineend = "butt", color = "grey70", linewidth = .1), axis.title.x = element_text(size = 13),
            axis.ticks.x = element_line(lineend = "butt", color = "grey70", linewidth = .1),
            legend.position = c(0.75,0.25), 
            legend.background = element_blank(), legend.key = element_blank(), legend.box.background = element_blank())
    
    #### Plot HYDE ####
    Plot.HYDE.Uz.1.5k <- Plot.HYDE(Extract = HYDE.Uz, Limits = c(-100, 1327), Log.scale = F, Facet_plots = F,
                                   Time.window = 100,  Leg.pos = c(0.6,0.2))
    
    #### Plot export ####
    Plot.compa <- Plot.compa + Plot.div + Plot.HYDE.Uz.1.5k + plot_layout(widths = c(10,2,2))
    H = 600; W = 1600; Save.plot ="Figures/Figure_4B.pdf" 
    ggsave(filename = Save.plot, Plot.compa, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
  }
  
  #### Figure 5 ####
  Fig5 = F
  if(Fig5 == T){
    Fazilman.biomiz$Age <- sort(MCWT.Faz.MP_ss_gf$MCWT$Age)
    Figure.5.plot <- CWM.FT.plot(MT = MCWT.Faz.MP_ss_gf,
                            Select.clim = c("MAAT", "MAP", "AI"), Mbiom = Fazilman.biomiz,
                            Select.model = c("MAT", "BRT"),
                            Select.trait = c("TRY_Height", "TRY_LeafArea", "TRY_SSD"),
                            Select.biom = c("COST", "WAST", "TUND", "CODE", "HODE", "CLMX", "COMX", "COCO"), Show.main.biomes = T, All.biomes = F,
                            Plot.x = "Age", Select.DB = NULL, Return.plot = T, Strat.plot = T, Dot.size = 1.3, Dot.alpha = .5, Pourc.biom = 50, Smooth.param = .24,
                            H = 800, W = 500, Display.legend = "bottom", Select.interv = 1000, Limites = c(0,10000), Lab.age.rotation = 25,
                            Zone.clim = c(50, 550, 650, 950, 1600, 1800, 1900, 2500, 3750, 4400, 5000, 7500, 8200, max(MCWT.Faz.MP_ss_gf$MCWT$Age)),
                            Temp.zone = c("C","W","C","W", "D", "Wt", "D"),
                            Name.zone = c("LIA", "WMP", "DACP", "RWP", "4.2k", "Mid-Hol.", "Early Hol."),
                            Zone.rotation = 30, Zone.clim.box = F, Panel.ecart = 0,
                            )
    
    Hist.diff <- MCWT.Faz.MP_ss_gf[[1]][c("TRY_SSD", "TRY_Height", "Age")]
    Hist.diff$TRY_Resist_arid <- (Hist.diff$TRY_SSD - Hist.diff$TRY_Height) / ((Hist.diff$TRY_SSD + Hist.diff$TRY_Height) / 2)
    Hist.d <- ggplot(Hist.diff, aes(x = Age, y = TRY_Resist_arid))+geom_col(width = 70, position = "dodgev", na.rm = T, fill = "grey30")
    Figure.5.plot <- Figure.5.plot[[1]] / Figure.5.plot[[2]] / Hist.d + plot_layout(heights = c(2, 40, 2))
    W = 500; H = 900; ggsave(filename = "Figures/Figure_5.pdf", Figure.5.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")
  }
  
  #### Figure 6 ####
  Fig6 = F
  if(Fig6 == T){
    #### Plot PCA GDGT ####
    GDGT.Fazilman.PCA <- GDGT.Fazilman[c(setdiff(grep("f.I", names(GDGT.Fazilman)), union(grep("7Me", names(GDGT.Fazilman)), grep("GDGT", names(GDGT.Fazilman)))), grep("Lith", names(GDGT.Fazilman)))]
    pca.Faz.GDGT <- PCA.XRF(GDGT.Fazilman.PCA, transp_OK = T, Scale.PCA = 3, return.plot = T, Legends.pos = "bottom",
                            Cluster.core = "Litho", 
                            Cluster.core.lab = "Fazilman lithological units", GDGT = T, Leg.nrows = 2,
                            Color.choice = Fazilman.clus.color, Show.lab.PCA = F,
                            Color.vectors = c("#87CA6A", "#75AADB", "#c67f05", "grey", "#E76D51", "#004266"),
                            Groupes = list("Alkalin conditions" = c("IIIa'", "IIa'"),
                                           "Warm conditions" = c("Ia", "Ib", "Ic", "IIb"),
                                           "Cold conditions" = c("IIb'", "IIc", "IIc'"),
                                           "Dry conditions" = c("IIb'", "IIIc'", "IIIb'"),
                                           "Wet conditions" = c("IIa", "IIIa")
                            ), 
                            Csv.sep = ",", Site.name = "(C)", Type.samples = "", Ellipse = T, Show.centroid = T,
                            # Save.plot = "Figures/Uzbekistan/GDGT/Fazilman/PCA_GDGT_Fazilman.pdf", H = 450, W = 570,
                            # Save.path = "Resultats/Uzbekistan/GDGT/Fazilman/PCA_GDGT_Fazilman.csv" 
                            )
    
    #### Plot scatterplots GDGT ####
    p1 <- Plot.BIT.IIIa(GDGT.Fazilman, #Manu.lim = c(0.2, 1.35, 0.9, 1), Vline = c(0.5,0.63), Hline = c(0.98, 0.96),
                        Cluster.core = "Litho", Color.choice = Fazilman.clus.color, Leg.pos = "none", 
                        Ellipse = T, Title = "(A)")
    
    p2 <- Plot.BIT.IIIa(GDGT.Fazilman, Select.param = c("Ib.Ia", "IIIa.IIa"), Leg.pos = "none", Title = "(B)",
                        Hline = c(0.7), #Manu.lim = c(0.6, 1.2, 0, 1.6),
                        Symbol.pos = c(.78,.78,.35), Ellipse = T,
                        # Symbol.path = "Figures/ACA/Trait/Workflow/Symboles/xml/Symbole_GDGT_2.xml",
                        Cluster.core = "Litho", Color.choice = Fazilman.clus.color)
    
    p3 <- Plot.BIT.IIIa(GDGT.Fazilman, Select.param = c("CI", "IIIa.IIa"),# Vline = c(0.5,0.63), Hline = c(0.35, 0.25),
                        Cluster.core = "Litho", Color.choice = Fazilman.clus.color, Leg.pos = "none", Title = "(C)")
    
    p4 <- Plot.BIT.IIIa(GDGT.Fazilman, Select.param = c("IIIa.IIa", "Ib.Ia"), Ellipse = T,
                        Cluster.core = "Litho", Color.choice = Fazilman.clus.color)
    
    p5 <- Plot.BIT.IIIa(GDGT.Fazilman, Select.param = c("MBTp5Me", "MBTp6Me"), Leg.pos = "none", Title = "(D)", Ellipse = F,
                        Manu.lim = c(0.15, 0.7, 0.15, 0.7), No.lines = T, Add.reg.lin = T, R2.pos = "topleft", Bisectrice = T,
                        Cluster.core = "Litho", Color.choice = Fazilman.clus.color)
    
    #### Plot delta(MBT5 / MBT6) ####
    A <- subset(GDGT.Fazilman, select = c(Age, MBTp5Me, MBTp6Me, Litho))
    A$Delta <- A$MBTp5Me - A$MBTp6Me
    A.mx <- subset(GDGT.Fazilman.MaxI, select = c(Age, MBTp5Me, MBTp6Me, Litho))
    A$Delta.max <- A.mx$MBTp5Me - A.mx$MBTp6Me
    A.mn <- subset(GDGT.Fazilman.MinI, select = c(Age, MBTp5Me, MBTp6Me, Litho))
    A$Delta.min <- A.mn$MBTp5Me - A.mn$MBTp6Me
    
    My_arrow <- data.frame(x1 = 2000, x2 = 2000, y1 = -0.05, y2 = min(A$Delta))
    # pDelta <- ggplot(A, aes(y = Delta, x = Age, color = Litho, group = 1))+
    pDelta <- ggplot(A, aes(y = Delta, x = Age, color = Litho, group = 1))+
      scale_y_reverse()+ ggtitle("(E)")+
      scale_x_continuous(breaks = seq(0, 10000, 1000))+
      scale_color_manual(values = setNames(Fazilman.clus.color,rev(unique(A$Litho))))+
      geom_ribbon(aes(ymin = Delta.min, ymax = Delta.max, x = Age), fill = "grey30", color = NA, alpha = .3)+
      geom_line(linewidth = .5, color = "grey30", orientation = "x")+
      geom_point(size = 2)+ xlab("Time (year cal BP)")+ ylab(expression(Delta(MBT*minute[5~Me]~";"~MBT*minute[6~Me])))+ 
      geom_vline(xintercept = 8200, linetype = "dashed", color = "grey30")+
      geom_vline(xintercept = 4200, linetype = "dashed", color = "grey30")+
      annotate("text", x = 2500, y = min(A$Delta)+0.05, label = "Wetter", color ="#6872B3", size = 6, hjust = 0, vjust = 1)+
      geom_segment(data = My_arrow, mapping = aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.5,"cm")), color = "#6872B3", linewidth = .8, lineend = "round")+
      theme(legend.position = "none", plot.background = element_blank(), panel.background = element_rect(NA, "black", linewidth = 1), panel.grid = element_blank(), plot.margin=unit(c(0,0,0,0),"cm"))
    
    #### Merge plots and export ####
    pca.CWT.full <- (p1 + p2)/(pca.Faz.GDGT + p5) / pDelta / guide_area() + plot_layout(guides = "collect")
    
    W = 600*0.9; H = 1150*0.9; Save.plot = "Figures/Figure_6.pdf"
    ggsave(filename = Save.plot, pca.CWT.full, width = W*0.026458333, height = H*0.026458333, units = "cm")  
  }
  
  #### Figure 7 ####
  Fig7 = F
  if(Fig7 == T){
    A <- MGDGT.change.name(List.GDGT = list(GDGT.Fazilman, GDGT.Fazilman.MinI, GDGT.Fazilman.MaxI),
                           Old.names = c("P.tetra", "P.penta", "P.hexa", "f.Ib", "f.IIa_5Me", "CI", "IR", "MAAT_Ensemble"),
                           New.names = c("Tetra.methyl", "Penta.methyl", "Hexa.methyl", "CBT.Ib.frac", "CBT.IIa.frac", "Rib.CI.community", "CBT.IR", "MAAT_BRT_Ensemble"))
    GDGT.Fazilman <- A[[1]]
    GDGT.Fazilman.MinI <- A[[2]]
    GDGT.Fazilman.MaxI <- A[[3]]
    
    Fazilman.plot.gdgt <- Plot.GDGT.clim(MGDGT = list(Fazilman = GDGT.Fazilman),
                                         MMin = list(Fazilman = GDGT.Fazilman.MinI),
                                         MMax = list(Fazilman = GDGT.Fazilman.MaxI),
                                         Surf.val = Cores.metadata,
                                         Cores = c("Fazilman"),
                                         Pclim = c("MBT", "methyl", "CBT",  "Rib", "MAAT", "AI"),
                                         Clim.lab = c("MBT indices", "Methylations", "Physicochemical changes", "Community changes", "MAAT-reconstructions (°C)", "AI"),
                                         Pourc.clim.param = c(1/10, 1/10, 2/10, 2/10, 2/10, 1/10),
                                         Mono.core = T, 
                                         Keep.models = c("MBTp6Me", "MBTp5Me", "Tetra.methyl", "Penta.methyl", "Hexa.methyl", "CBT.Ib.frac", "CBT.IIa.frac", "Rib.CI.community", "CBT.IR",
                                                         "CBTp", "IR6Me", "AI_BRT_ACA_CB_LS", "AI_Ensemble", "AI_BRT_ACA_LS", "Rib", "MAAT_Chen_Tjk", "MAAT_DJ_5Me", "MAAT_mr_Yang",
                                                         "MAAT_LSun", "MAAT_mr_DJ", "MAAT_mrs_DJ", "MAAT_BRT_ACA_LS", "MAAT_BRT_Ensemble", "MAAT_BRT_ACA_CB_LS"),
                                         Label.group = c("community", "frac", "Crenar", "MAAT_BRT"),
                                         Keep.other.grp = T,
                                         Age.select = "Age", Surf.val.max.age = 1000,
                                         Pollen.plot.merge = F, Manual.vlines = 8200, Repel.repoussage = 10,
                                         Annotation = F, Show.x.axis = T, Shape.group = c(20,20,20,20,20,20,20),
                                         Facet.T = T, Show.repels = T, Repel.x = 12700, Lab.clim.angle = 0,
                                         Facet.scale.egal = F, Display.legends = "none", Panel.box = T,
                                         Select.interv = 1000, 
                                         Zone.clim = c(50, 550, 650, 950, 1600, 1800, 1900, 2500, 3300, 3900, 4000, 4400, 5000, 7650, 7650, 8950, 9500 , max(GDGT.Fazilman2$Age)),  #Feng et al., 2006
                                         Temp.zone = c("C","W","C","W","C", "W", "C", "C", "G"), Zone.clim.dashed = F,
                                         Name.zone = c("LIA", "WMP", "DACP", "RWP", "3.5ky", "4.2ky", "Lak1",  "Lak2", "LS"),
                                         Save.plot = "Figures/Figure_7.pdf",
                                         H = 1900, W = 680)  
  }
  
  #### Figure 8 ####
  Fig8 = F
  if(Fig8 == T){
    #### Plot Obs vs. Pred (buffer) ####
    M.MAAT.GDGT <- Msurf.Faz[,grep("^MAAT", colnames(Msurf.Faz))]
    p <- Residual.plot(MML = NULL, Msurf = M.MAAT.GDGT[c(2,8,9,14,20)], Dot.size = 2.5, Res.right = F, Add.RMSE = T, Nrow.plot = 1, Add.nb = F, 
                       Residual.lims = c(-5, 10), Annot.position = "topleft", CV.method = "self",
                       Mclim = Mclim.Faz, Annot = paste("D", seq(1:5), sep = ''), Hexa = F)
    
    #### Récupération des points proches de Fazilman ####
      M1 <- Compar.Surf.Calib(Msurf.Faz, Mclim.Faz, 
                              Select.clim = "MAAT",
                              Clim.display = "Anomaly", 
                              H.arrow = F, 
                              Calib.to.keep = c("MAAT_BRT_ACA", "MAAT_Chen_Tjk", "MAAT_mr_DJ", "MAAT_DJ_5Me", "MAAT_BRT_ACA_LS"),
                              Bold_dots = "MAAT_BRT_ACA_LS", Bold.width = 2,
                              Barres = c(3.7), Arrow.width = 0.39, Add.frame = T,
                              Order.by = Mclim.Faz["Dist"], Lab.angle = 35, Lab.size = 8,
                              Box.order.mean = T, Leg.row = 1, Break.line.annot = F,
                              Lab.area = c("< 10 km", "> 10 km", "> 100 km"),
                              Pt.size = .5)
    
      #### Import data SIG ####
        Eurasia_map <- c("Mongolia", "Russia", "Spain", "France", "Italy", "Greece", "Germany", "Finland", "Bhutan", "Bangladesh",
                         "Czech Republic", "Denmark", "Kosovo", "Vietnam", "Laos", "Japan", "Nepal", "India", "Myanmar",
                         "Latvia","Lithuania", "Estonia", "Belarus", "Romania", "Bulgaria", "Hungary", "Austria", "Croatia",
                         "Albania", "Serbia", "Slovenia", "Slovakia", "Bosnia", "Montenegro", "UK", "Ireland", "Moldova", "Macedonia",
                         "Norway", "Sweden", "Turkey", "Belgium", "Uzbekistan", "Tajikistan", "Syria", "Israel", "Jordan", "Pakistan",
                         "Kazakhstan","Turkmenistan", "Ukraine", "Poland", "Portugal", "Switzerland", "Kyrgyzstan", "Morocco",
                         "China", "Iran", "Armenia", "Georgia", "Afghanistan", "Iraq", "Azerbaijan")
        Eurasia_map <- map_data("world", region = Eurasia_map)
      
      #### Plot k-mean map ####
      pmap <- ggplot(Mclim.Faz, aes(y = Latitude, x = Longitude))+
        new_scale_color()+ new_scale_fill()+
        geom_polygon(data = Eurasia_map, aes(x=long, y=lat, group = group), alpha = 1, fill = NA, color = "grey10", size = 0.3)+
        geom_point(data = Meco.Faz[row.names(Meco.Faz) == "MUZT6C02",], size = 8, fill = "darkorange", shape = 21)+
        geom_point(size = 4, alpha = 1, aes(fill = Dist), shape = 21)+
        scale_fill_continuous(name = "Distance to Fazilman (km)")+
        ggforce::geom_ellipse(data = Meco.Faz[row.names(Meco.Faz) == "MUZT6C02",], 
                              aes(y0 = Meco.Faz[row.names(Meco.Faz) == "MUZT6C02",]$Latitude, x0 = Meco.Faz[row.names(Meco.Faz) == "MUZT6C02",]$Longitude, a = 2.6, b = 2.3, angle = 0),
                              linetype = 2, color = "grey20")+
        coord_fixed(xlim = c(62, 71), ylim = c(37, 43)) +
        labs(x = c("Longitude (°)"), y = c("Latitude (°)"))+ 
        theme(legend.position = "top", panel.background = element_rect(fill = NA, colour = "black"), panel.grid = element_blank(),
              axis.title = element_text(margin = ggplot2::margin(t = -2, r = 0, b = 0, l = 0), size = 15), legend.background = element_blank())
    
    #### Full plot ####
      p <- (pmap + M1) / p
      pp = 0.95
      H = 950*pp; W = 1500*pp; Save.plot ="Figures/Figure_8.pdf"
      ggsave(filename = Save.plot, p, width = W*0.026458333, height = H*0.026458333, units = "cm")
  }
  
  #### Figure 9 ####
  Fig9 = F
  if(Fig9 == T){
    Figure.9.plot <- Plot.FT.summary(FT = c(Fazilman = c(Fazilman.TUSDB, Fazilman.ACADB)),
                                   RMSE.barres = c(TUSDB = c(MAT = MAT.TUSDB, WAPLS = WAPLS.TUSDB, BRT = BRT.TUSDB), 
                                                   ACADB = c(MAT = MAT.ACADB, WAPLS = WAPLS.ACADB, BRT = BRT.ACADB)),
                                   Cores = c("Fazilman"),
                                   Pclim = c("MAAT", "MPCOQ", "AI"),
                                   Label.group = c("TUSDB", "COSTDB", "ACADB"),
                                   Surf.val = Cores.metadata,
                                   Anomaly = F, Mono.core = T, Merge.legends = T, Legend.pos = "bottom",
                                   GDGT.plot.merge = F, Add.lim.space = F, Smooth.T = T, Condensed = T,
                                   Clim.lab = c("(A) MAAT (°C)", expression("("*B*")"~MPCOQ~(mm.yr^-1)), "(C) AI", "MTWAQ (°C)"),
                                   Model = c("WAPLS", "MAT", "BRT"), 
                                   Select.interv = 1000, 
                                   Limites = c(-100, 10000),
                                   Manual.y.val = list(MAAT.Fazilman = c(-2, 0, 2, 4, 6, 8, 10), MPCOQ.Fazilman = c(0, 50, 100, 150, 200, 250, 300), AI.Fazilman = c(2000, 3000, 4000, 5000)),
                                   Zone.clim = c(50, 550, 650, 950, 1600, 1800, 1900, 2500, 3750, 4400, 5000, 7500, 8200, max(Fazilman.TUSDB$WAPLS.TUSDB$Age)),
                                   Temp.zone = c("C","W","C","W", "D", "Wt", "D"), Cores.lab = "Pollen-based climate reconstructions",
                                   Name.zone = c("LIA", "WMP", "DACP", "RWP", "4.2k", "Mid-Hol.", "Early Hol."), 
                                   Temp.col.alpha = 0.1, Zone.clim.box = F, Name.zone.angle = 0, 
                                   Save.plot = "Figures/Figure_9.pdf",
                                   H = 800, W = 1200)
  }
  #### Figure 10 ####
  Fig10 = F
  if(Fig10 == T){
    #### Pollen-based Ensemble models ####
    FT.Ens <- Pollen.Ensemble.RMSE(MPal = c(Fazilman.TUSDB, Fazilman.ACADB),
                                   Msurf = list("BRT.TUSDB"   = BRT.TUSDB,  "BRT.ACADB"  = BRT.ACADB, 
                                                "MAT.TUSDB"   = MAT.TUSDB,  "MAT.ACADB"  = MAT.ACADB, 
                                                "WAPLS.TUSDB" = WAPLS.TUSDB,"WAPLS.ACADB"= WAPLS.ACADB),
                                   Limites = c(-100, 10000),
                                   Remove.calib = c("MPCOQ.MAT.TUSDB", "AI.MAT.TUSDB", "MAAT.MAT.TUSDB",
                                                    "AI.MAT.ACADB", "AI.WAPLS.ACADB", "MAAT.BRT.ACADB"),
                                   Surf.val = Cores.metadata, Surf.val.max.age = 1000,
                                   Clim.lab = c("(A) MAAT (°C)", expression("("*B*")"~MPCOQ~(mm.yr^-1)), "(C) AI"),
                                   Manual.vlines = c(50, 4200, 5820, 6300, 8200),
                                   Select.interv = 500, Lab.age.size = 10, Temp.col.alpha = 0.2, Dot.size = 2.7,
                                   Manual.y.val = list(MAAT.Fazilman = c(2, 4, 6, 8, 10, 12), MPCOQ.Fazilman = c(50,75,100,125,150), AI.Fazilman = c(2000,2500,3000,3500,4000)),
                                   GDGT.plot.merge = T, Smooth.T = T, Smooth.sd = F, Smooth.param = 0.25,
                                   Name.core = "Fazilman", return.plot = T, Show.box = F,
                                   Cores.lab = c("Pollen-based \n climate reconstructions"),
                                   Plot.x = "Age", Keep.clim = c("MAAT", "MPCOQ", "AI"))
    
    #### brGDGT-based Ensemble models ####
    Faz10k <- Plot.GDGT.clim(MGDGT = list(Fazilman = GDGT.Fazilman),
                             MMin = list(Fazilman = GDGT.Fazilman.MinI),
                             MMax = list(Fazilman = GDGT.Fazilman.MaxI),
                             Surf.val = Cores.metadata, Age.select = "Age", Cores.lab = c("brGDGT-based \n climate reconstructions"),
                             Surf.val.max.age = 1000,
                             Pclim = c("MAAT", "MPCOQ", "AI"),
                             Clim.lab = c("MAAT (°C)", expression(paste(MPCOQ ~~ (mm.yr^1))), "AI"), 
                             Title.age = "Time (year cal BP)", 
                             Smooth.param = 0.25, Smooth.SD = F,
                             Label.group = c("Ens", "BRT"), Shape.group = c(1,1,2,5),
                             Keep.models = c("MPCOQ_BRT_ACA_LS", "MAF_Ensemble", "AI_Ensemble", "MAAT_Ensemble"),
                             Manual.y.val = list(MAAT.Fazilman = c(2,4,6,8,10,12), MPCOQ.Fazilman = c(50, 75, 100, 125, 150), AI.Fazilman = c(1500, 2000, 2500, 3000, 3500, 4000)),
                             Mono.core = F, Multi.clim = T, Add.lim.space = F,  Pollen.plot.merge = T, Annotation = F, 
                             Keep.other.grp = F, Anomaly = F, Facet.T = F, Facet.size.egal = T, Facet.scale.egal = F, Select.interv = 500,
                             Limites = c(-100, 10000), Show.repels = F, Manual.vlines = c(50, 4200, 5820, 6300, 8200),
                             Temp.zone = c("W", "D","Wt","D","Wt","C", "G"), #, "W"
                             Panel.box = F, Color.by.proxy = T, Show.x.axis = T,
                             Temp.col.alpha = 0.1, Dot.size = 2.7, Display.legends = "none",
                             Save.plot = "Figures/Figure_10.pdf",
                             H = 400, W = 1500)
    
    #### Chronological diagram ####
    Chrono.MAP <- Timeline.build(
      Datatime = Datatime.MAP, Limites = c(-100, 10000), Leg.pos = "bottom",
      Add.AC.BC = T, Label.size = 4.2, Select.interv = 1000, Axis.pos.y = 0 )
    
    Chrono.MAAT <- Timeline.build(
      Datatime = Datatime.MAAT, Limites = c(-100, 10000), Leg.pos = "bottom",
      Add.AC.BC = T, Label.size = 4.2, Select.interv = 1000, Axis.pos.y = 0)
    
    #### Full plot ####
    Save.plot = "Figures/Figure_10.pdf"
    
    layout <- 
      "
            AAAAAAAAAADDDDDDDDDDGGGGGGGGG
            AAAAAAAAAADDDDDDDDDDGGGGGGGGG
            BBBBBBBBBBEEEEEEEEEEHHHHHHHHH
            BBBBBBBBBBEEEEEEEEEEHHHHHHHHH
            BBBBBBBBBBEEEEEEEEEEHHHHHHHHH
            CCCCCCCCCCFFFFFFFFFFIIIIIIIII
            CCCCCCCCCCFFFFFFFFFFIIIIIIIII
            CCCCCCCCCCFFFFFFFFFFIIIIIIIII
             "
    pFaz10k <- Chrono.MAAT + FT.Ens[[1]] + Faz10k[[1]] + Chrono.MAP + FT.Ens[[2]] + Faz10k[[2]] + plot_spacer() + FT.Ens[[3]] + Faz10k[[3]] + 
      plot_layout(design = layout, guides = "collect") & theme(plot.margin = unit(c(0,0,0,0),"cm"), legend.position = "bottom")
    
    W = 1800; H = 800
    ggsave(filename = Save.plot, pFaz10k, width = W*0.026458333, height = H*0.026458333, units = "cm")
  }
  
  #### Figure S4 ####
  FigS4 = F
  if(FigS4 == T){
    Figure.S5.plot <- DP_plot_algue(Fazilman.AlF, MP.Fazilman.Age,
                                    CONISS = F, Time.window = 500, Nzone = 5,
                                    Save.plot = "Figures/Figure_S4.pdf",
                                    H = 1100, W = 1900)
  }
  
  #### Figure S5 ####
  FigS5 = F
  if(FigS5 == T){
    Figure.S5.plot <- DP_pol(Fazilman.PoF, MP.Fazilman.Age, Plot.x = "Age",
                             Mpol_conc = Fazilman.PoIx,
                             Path.index = "Import/Dictionary/Pollen_table.csv",
                             Seuil.Pour = 0.0001, Time.window = 500, nlake = NULL, 
                             CONISS = F, Nzone = 10, AP.NAP = F, Only.AP = T, Sort.taxon = "Auto.mean", 
                             Pollen.diversity = T, Spore.diversity = F, Auto.Pollen.lab = T, 
                             Pol.influx.total = T, Spore.influx = F, Cortege = F, Richesse.type = F,
                             Save.plot = "Figures/Figure_S5.pdf",
                             H = 900, W = 2800
    )
  }
  
  #### Figure S6 ####
  FigS6 = F
  if(FigS6 == T){
    Fazilman.biomiz$Age <- sort(MCWT.Faz.MP_ss_gf$MCWT$Age)
    Figure.S6.plot <- CWM.FT.plot(Mbiom = Fazilman.biomiz,
                            Show.main.biomes = T, Plot.x = "Age", Select.DB = NULL, Return.plot = T, Strat.plot = T, Dot.size = 1.3, Dot.alpha = .5,
                            Pourc.biom = 40, Smooth.param = .24, Biom.groups = T, All.biomes = T,
                            H = 1200, W = 500, Display.legend = "bottom", Select.interv = 1000, Limites = c(0,10000),
                            Save.plot = "Figures/Figure_S6.pdf")}
  #### Figure S7 ####
  FigS7 = F
  if(FigS7 == T){
    Figure.S7.plot <- CWM.FT.plot(MT = MCWT.Faz.MP_ss_gf,
                            Select.trait = c("TRY_Height", "TRY_LeafArea", "TRY_SeedMass", "TRY_SLA", "TRY_SSD", "TRY_LeafN"),
                            Plot.x = "Age", Return.plot = T, Strat.plot = T, Dot.size = 1.3, Dot.alpha = .5, Smooth.param = .24,
                            H = 1000, W = 500, Display.legend = "none", Select.interv = 1000, Limites = c(0,10000),
                            Save.plot = "Figures/Figure_S7.pdf")}
  
  #### Figure S8 ####
  FigS8 = F
  if(FigS8 == T){
    #### Test et vérifications - Fazilman #### 
    test.randomTF = F
    if(test.randomTF == T){
      Faz.rTF.TUSDB <- Plot.randomTF(MPsurf = TUSDB.MP, MPpaleo = MP_Fazilman.conv,  Mclim = Uz.clim[randomTF.clim], Database = "TUSDB", Lake = "Fazilman", Plot.MAT = T, Plot.WAPLS = T, 
                                     Save.path = "Resultats/Uzbekistan/Pollen/Func_trans/Fazilman/Fazilman_randomFT_TUSDB.Rds", return.plot = T, Plot.RF = F, Plot.BRT = T,
                                     H = 500, W = 2000, Save.plot = "Figures/Uzbekistan/Pollen/Func_trans/Fazilman/randomTF/Fazilman_randomTF_TUSDB.pdf")
      saveRDS(Faz.rTF.TUSDB, "Resultats/Uzbekistan/Pollen/Func_trans/Fazilman/Fazilman_randomFT_TUSDB_plot.Rds")
    }
    
  }
}