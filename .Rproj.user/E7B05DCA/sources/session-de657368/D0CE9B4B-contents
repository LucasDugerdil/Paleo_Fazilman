#### Functions ####
Diversite_pollen <- function(Mpol){
  MDivers <- data.frame(Diversite=0)
  DivNul <- colSums(Mpol, na.rm = F)
  DivNul[is.na(DivNul)] <- 0
  for (nsite in 1:length(Mpol)){
    counter = 0
    if(DivNul[nsite] == 0){MDivers[nsite]=0}
    else{
      for (ntaxon in 1:length(Mpol[,1])){
        if(Mpol[ntaxon,nsite] != 0){
          counter = counter + 1
          MDivers[nsite]=counter
        }}}}
  names(MDivers) <- names(Mpol)
  return(MDivers)
}

CWT.calculation <- function(MT, MP, Mclim = NULL, MPS.ACA.Biom = NULL, Accep.seuil, Remove.biom, Add.CWV = F){
  #### Settings ####
  if(missing(Remove.biom)){Remove.biom = NULL}
  
  #### Clean matrix ####
  if(any(colSums(MP) < 95)){
    print("Recalculate the FA for the species matrix.")
    Keep.real.names <- row.names(MP)
    MP <- data.frame(t(MP))
    MP <- MP/rowSums(MP)*100
    MP <- data.frame(t(MP))
    row.names(MP) <- Keep.real.names
  }
  
  MP$id <- row.names(MP)
  names(MT) <- gsub("X", "", names(MT))
  if(any(grepl("TRY", names(MT)) == T) == F){names(MT) <- paste("TRY", names(MT), sep = "_")}
  names(MT)[1] <- "id"
  MT <- MT[which(MT$id %in% MP$id),]
  # /!\ missing taxa ! 
  Missing.taxa <- setdiff(MP$id, MT$id)
  if(length(Missing.taxa) > 0){
    print("The following taxa are missing from the trait matrix:")
    print(Missing.taxa)
  }
  
  Missing.raw <- MT[setdiff(MP$id, MT$id),]
  Missing.raw$id <- setdiff(MP$id, MT$id)
  MT <- rbind(MT, Missing.raw)
  MT <- as_tibble(lapply(MT, function(x){x[is.nan(x)] <- NA ; x}))
  
  #### Merge MP + MT ####
  MPT <- merge(MT, MP, by = "id", all = T)
  MPT <- melt(MPT, id = names(MT))
  names(MPT)[names(MPT) == "variable"] <- "Site"
  names(MPT)[names(MPT) == "value"] <- "FA"
  
  #### Calculation of the CWT ####
  MCWT <- data.frame(Site = (unique(MPT["Site"])))
  MCWT.stat <- setNames(data.frame(NA,NA,NA), c("Trait", "Pour.sites.kept", "N.site"))
  MCWT.stat <- MCWT.stat[-1,]
  Tot.site.nb <- length(levels(MPT$Site))
  
  if(Add.CWV == T){MCWV <- MCWT; MCWV.stat = MCWT.stat}
  
  #### Main loop ####
  for(i in grep("TRY", names(MPT))){
    Trait.treat.i <- names(MPT)[i]
    A <- MPT[c("id", Trait.treat.i, "Site", "FA")]
    A$FA[is.na(A[[Trait.treat.i]])] <- NA
    if(all(is.na(A$FA)) == T){next}
    all_abund = aggregate(FA ~ Site, A, sum)
    colnames(all_abund)[2] = "tot_abund"
    Keep.sites <- all_abund[all_abund[2] > Accep.seuil,]
    
    #### Calculation CWM ####
    if(nrow(Keep.sites) > 0){
      N.site <- length(Keep.sites$tot_abund)
      Pourc.site.up.seuil <- round(length(Keep.sites$tot_abund)/Tot.site.nb, digits = 2)
      MCWT.stat[i,] <- c(Trait.treat.i, Pourc.site.up.seuil, N.site)
      A <- A[which(A$Site %in% Keep.sites$Site),]
      A <- merge(A, all_abund, by = "Site")
      # A$FA <- scale(A$FA/A$tot_abund)
      A$FA <- A$FA/A$tot_abund
      XX <- aggregate(FA * eval(parse(text = Trait.treat.i)) ~ Site, A, sum, na.rm = T)
    }
    else{
      XX <- data.frame(Site = NA, X = "")
    }
    names(XX)[2] <- Trait.treat.i
    MCWT <- left_join(MCWT, XX, by = "Site")
    
    
    #### Add CWV ####
    if(Add.CWV == T){
      if(nrow(Keep.sites) > 0){
        
        A.var <- full_join(A, setNames(XX, c("Site", "CWM")), by = "Site")
        A.var$T.var <- (A.var[[names(A.var)[grep("TRY_", names(A.var))]]] - A.var$CWM)^2
        XX.var <- aggregate(FA * T.var ~ Site, A.var, sum, na.rm = T)
      }
      else{
        XX.var <- data.frame(Site = NA, X = "")
      }
      names(XX.var)[2] <- Trait.treat.i
      MCWV <- left_join(MCWV, XX.var, by = "Site")
    }
    
  }
  
  MCWT.stat <- MCWT.stat[-1,]
  MCWT.stat[nrow(MCWT.stat)+1,] <- c("Average", round(mean(as.numeric(MCWT.stat$Pour.sites.kept)), digits = 2), round(mean(as.numeric(MCWT.stat$N.site)), digits = 0))
  
  #### Merge CWT + climat ####
  if(is.null(Mclim) == F){
    Mclim[["Site"]] <- row.names(Mclim)
    MCWT.clim = merge(MCWT, Mclim, by = "Site")    # On fusionne les matrices CWM et CLIMAT
    if(Add.CWV == T){MCWV.clim = merge(MCWV, Mclim, by = "Site")}
    
    if(is.null(MPS.ACA.Biom) == F){
      MPS.ACA.Biom[["Site"]] <- row.names(MPS.ACA.Biom)
      MCWT.clim = merge(MCWT.clim, MPS.ACA.Biom, by = intersect(names(MCWT.clim), names(MPS.ACA.Biom)))    # On fusionne les matrices CWM et CLIMAT
      if(Add.CWV == T){MCWV.clim = merge(MCWV.clim, MPS.ACA.Biom, by = intersect(names(MCWV.clim), names(MPS.ACA.Biom)))}
    }
    if(is.null(Remove.biom) == F){
      MCWT.clim <- MCWT.clim[setdiff(seq(1,nrow(MCWT.clim)), which(MCWT.clim$Biome %in% Remove.biom)),]
      if(Add.CWV == T){MCWV.clim <- MCWV.clim[setdiff(seq(1,nrow(MCWV.clim)), which(MCWV.clim$Biome %in% Remove.biom)),]}
    }
  }
  
  #### Export ####
  if(Add.CWV == F){
    if(is.null(Mclim) == F){
      Lexport <- list(MCWT = MCWT.clim, MCWT.stat = MCWT.stat, Missing.taxon = Missing.taxa)
      if(is.null(Mclim) == T){
        Lexport <- list(MCWT = MCWT, MCWT.stat = MCWT.stat, Missing.taxon = Missing.taxa)}
    }}
  else{
    if(is.null(Mclim) == F){
      Lexport <- list(MCWT = MCWT.clim, MCWV = MCWV.clim, MCWT.stat = MCWT.stat, Missing.taxon = Missing.taxa)
      if(is.null(Mclim) == T){
        Lexport <- list(MCWT = MCWT, MCWV = MCWV, MCWT.stat = MCWT.stat, Missing.taxon = Missing.taxa)}
      
    }}
  return(Lexport)
}

Biomization_cal <- function(Tax2PFT = NULL, PFT2biom = NULL, Pollen = NULL, Display.warning = T, bioindic = F, Verbose = T, 
                            Posterior.PFT.correction = T, PDF.reassign = NULL, Save.path = NULL) {
  #### My algo ####
  if(bioindic == F){
    #### Table PFT cleaning ####
    Tax2PFT[Tax2PFT == 0] <- NA
    Remove.pft <- Tax2PFT$RECN[!Tax2PFT$RECN %in% names(Pollen)]
    Tax2PFT <- Tax2PFT[!Tax2PFT$RECN %in% Remove.pft,]
    
    #### Table Pollen Cleaning ####
    if(Display.warning == T){
      Missing.tax <- names(Pollen)[!names(Pollen) %in% Tax2PFT$RECN]
      print("The following pollen-types are missing from the Taxa-to-PFT table & have been removed from the pollen matrix.")  
      print(Missing.tax)  
    }
    
    Pollen <- Pollen[, names(Pollen) %in% Tax2PFT$RECN]
    Pollen <- Pollen / rowSums(Pollen)*100
    Pollen <- Pollen[, match(Tax2PFT$RECN,names(Pollen))]
    
    #### Build result matrix ####
    Biom.output <- setNames(data.frame(matrix(0, nrow = nrow(Pollen), ncol = 1 + length(PFT2biom$RECN), byrow = T)), c("Sites", PFT2biom$RECN))
    Biom.output[,1] <- row.names(Pollen)
    
    PFT2biom <- cbind(row.names(PFT2biom), PFT2biom)
    PFT.output <- setNames(data.frame(matrix(0, nrow = nrow(Pollen), ncol = ncol(PFT2biom)-1, byrow = T)), c("Sites", names(PFT2biom)[-c(1,2)]))
    PFT.output[,1] <- row.names(Pollen)
    List.PFT <- names(PFT2biom)[-c(1,2)]
    
    #### PFT score for each sites ####
    if(Verbose == T){print("PFT scores calculation."); pb <- txtProgressBar(min = 0, max = nrow(Pollen)*length(List.PFT), style = 3, label = "PFT scores."); k<-0}
    for(i in 1:nrow(Pollen)){
      for(j in 1:length(List.PFT)){
        if(Verbose == T){k<-k+1; setTxtProgressBar(pb, k)}
        PFT.j <- List.PFT[j]
        PFT.j.val <- unlist(Tax2PFT[names(Tax2PFT) == PFT.j], use.names = F)
        Pollen.i <- unlist(Pollen[i,], use.names = F)
        X <- Pollen.i[Pollen.i >= PFT.j.val]
        X <- sum(sqrt(X[!is.na(X)]), na.rm = T)
        PFT.output[i,j+1] <- round(X, digits = 1)
      }
    }
    if(Verbose == T){close(pb)}
    
    #### Biomization calculation ####
    if(Verbose == T){print("Biome scores calculation."); pb <- txtProgressBar(min = 0, max = nrow(Pollen)*length(PFT2biom$RECN), style = 3); k <- 0}
    for(i in 1:nrow(Pollen)){
      for(j in 1:length(PFT2biom$RECN)){
        if(Verbose == T){k<-k+1; setTxtProgressBar(pb, k)}
        biome.j <- PFT2biom$RECN[j]
        biome.j <- unlist(PFT2biom[PFT2biom$RECN == biome.j,], use.names = F)
        biome.j[biome.j == "-"] <- NA
        biome.j[biome.j == "+"] <- NA
        biome.j[c(1,2)] <- NA
        PFT.j <- round(sum(PFT.output[i,which(biome.j == toupper(biome.j))-1], na.rm = T), digits = 1)
        Biom.output[i,j+1] <- PFT.j
      }
    }
    if(Verbose == T){close(pb)}
    
    #### Select best biome score by site ####
    Biom.output <- cbind(Biompol = names(Biom.output)[apply(Biom.output, 1, which.max)], Biom.output)
    PFT.output <- cbind(Biompol = names(Biom.output)[apply(Biom.output, 1, which.max)], PFT.output)
    
    #### Correction PFTs composite MAJUSCULE ####
    PFT.precor.t = F
    if(PFT.precor.t == T){
      PFT2biom2 <- data.frame(apply(PFT2biom, 2, function(x) gsub("-", NA, x)))
      PFT2biom2 <- data.frame(apply(PFT2biom2, 2, function(x) gsub("\\+", NA, x)))
      PFT2biom2 <- data.frame(apply(PFT2biom2, 2, function(x) ifelse(x == toupper(x), x, NA)))
      
      for(i in 1:nrow(Pollen)){
        print(paste("Biomes:", PFT.output$Biompol[i]))
        biome.i <- PFT2biom2[which(PFT2biom2$RECN == PFT.output$Biompol[i]),]
        
        for(j in 1:length(List.PFT)){
          if(is.na(biome.i[[j+2]]) == F & is.na(List.PFT[j]) == F){
            PFT.to.remove <- biome.i[[j+2]]
            PFT.output[i, names(PFT.output) == PFT.to.remove] <- PFT.output[i,j+2] + PFT.output[i, names(PFT.output) == PFT.to.remove]
            PFT.output[i,j+2] <- 0
          }}}
      
      Biom.output <- Biom.output[-c(1)]
      PFT.output <- PFT.output[-c(1)]
      for(i in 1:nrow(Pollen)){
        for(j in 1:length(PFT2biom$RECN)){
          biome.j <- PFT2biom$RECN[j]
          biome.j <- unlist(PFT2biom[PFT2biom$RECN == biome.j,], use.names = F)
          biome.j[biome.j == "-"] <- NA
          biome.j[biome.j == "+"] <- NA
          biome.j[c(1,2)] <- NA
          PFT.j <- sum(PFT.output[i,which(biome.j == toupper(biome.j))-1], na.rm = T)
          Biom.output[i,j+1] <- PFT.j
        }}
      
      Biom.output <- cbind(Biompol = names(Biom.output)[apply(Biom.output, 1, which.max)], Biom.output)
      PFT.output <- cbind(Biompol = names(Biom.output)[apply(Biom.output, 1, which.max)], PFT.output)
    }
    
    #### Securité + / - ####
    Security.t = T
    if(Security.t == T){
      if(Verbose == T){print("Cold/warm security"); pb <- txtProgressBar(min = 0, max = nrow(Pollen), style = 3); k <- 0}
      
      for(i in 1:nrow(Pollen)){
        if(Verbose == T){k<-k+1; setTxtProgressBar(pb, k)}
        
        if(PFT.output$Biompol[i] == "TUND"){
          Shift.PFT <- names(PFT2biom)[which(PFT2biom[PFT2biom$RECN == "TUND",] == "+")]
          if(sum(PFT.output[i,which(names(PFT.output) %in% Shift.PFT)]) > 0){
            PFT.output$Biompol[i] <- "COST"
            Biom.output$Biompol[i] <- "COST"}}
        
        if(PFT.output$Biompol[i] == "COST"){
          Shift.PFT <- names(PFT2biom)[which(PFT2biom[PFT2biom$RECN == "COST",] == "+")]
          if(sum(PFT.output[i,which(names(PFT.output) %in% Shift.PFT)]) > 0){
            PFT.output$Biompol[i] <- "WAST"
            Biom.output$Biompol[i] <- "WAST"}}
        
        if(PFT.output$Biompol[i] == "WAST"){
          Shift.PFT <- names(PFT2biom)[which(PFT2biom[PFT2biom$RECN == "WAST",] == "-")]
          if(sum(PFT.output[i,which(names(PFT.output) %in% Shift.PFT)]) > 0){
            PFT.output$Biompol[i] <- "COST"
            Biom.output$Biompol[i] <- "COST"}}
        
        if(PFT.output$Biompol[i] == "CODE"){
          Shift.PFT <- names(PFT2biom)[which(PFT2biom[PFT2biom$RECN == "CODE",] == "+")]
          if(sum(PFT.output[i,which(names(PFT.output) %in% Shift.PFT)]) > 0){
            PFT.output$Biompol[i] <- "HODE"
            Biom.output$Biompol[i] <- "HODE"}}
        
        if(PFT.output$Biompol[i] == "HODE"){
          Shift.PFT <- names(PFT2biom)[which(PFT2biom[PFT2biom$RECN == "HODE",] == "-")]
          if(sum(PFT.output[i,which(names(PFT.output) %in% Shift.PFT)]) > 0){
            PFT.output$Biompol[i] <- "CODE"
            Biom.output$Biompol[i] <- "CODE"}}
      }
      if(Verbose == T){close(pb)}
    }
    
    #### Correction PFTs composite minuscule ####
    if(Posterior.PFT.correction == T){
      PFT2biom2 <- data.frame(apply(PFT2biom, 2, function(x) gsub("-", NA, x)))
      PFT2biom2 <- data.frame(apply(PFT2biom2, 2, function(x) gsub("\\+", NA, x)))
      PFT2biom2 <- data.frame(apply(PFT2biom2, 2, function(x) toupper(x)))
      
      for(i in 1:nrow(Pollen)){
        biome.i <- PFT2biom2[which(PFT2biom2$RECN == PFT.output$Biompol[i]),]
        for(j in 1:length(List.PFT)){
          if(is.na(biome.i[[j+2]]) == F & is.na(List.PFT[j]) == F){
            PFT.to.remove <- biome.i[[j+2]]
            PFT.output[i, names(PFT.output) == PFT.to.remove] <- round(PFT.output[i,j+2] + PFT.output[i, names(PFT.output) == PFT.to.remove], digits = 3)
            PFT.output[i,j+2] <- 0
          }}}}
    
    #### Select best biome score by site ####
    
    #### Export results ####
    row.names(PFT.output) <- PFT.output$Sites
    row.names(Biom.output) <- Biom.output$Sites
    Biom.output <- Biom.output[-c(2)]
    PFT.output <- PFT.output[-c(2)]
    names(PFT.output)[names(PFT.output) == "Biompol"] <- "BIOMPOL"
    names(Biom.output)[names(Biom.output) == "Biompol"] <- "BIOMPOL"
    
  }
  #### Bioindic (Guiot) ####
  if(bioindic == T){
    library("bioindic")
    #### Table PFT cleaning ####
    Remove.pft <- row.names(Tax2PFT)[!row.names(Tax2PFT) %in% names(Pollen)]
    Tax2PFT <- Tax2PFT[!row.names(Tax2PFT) %in% Remove.pft,]
    
    taxcar <- Tax2PFT$Temp
    Tax2PFT <- subset(Tax2PFT, select = -c(Temp))
    Tax2PFT[is.na(Tax2PFT)] <- 0
    
    
    
    #### Table Pollen Cleaning ####
    if(Display.warning == T){
      Missing.tax <- names(Pollen)[!names(Pollen) %in% row.names(Tax2PFT)]
      print("The following pollen-types are missing from the Taxa-to-PFT table & have been removed from the pollen matrix.")  
      print(Missing.tax)  
    }
    
    Pollen <- Pollen[, names(Pollen) %in% row.names(Tax2PFT)]
    Pollen <- Pollen / rowSums(Pollen)*100
    Pollen <- Pollen[, match(row.names(Tax2PFT), names(Pollen))]
    
    #### Calcul biome ####
    pft <- pftscores(ass = Pollen, 
                     tax2pft = Tax2PFT, 
                     taxcar = taxcar, 
                     pftassign = PDF.reassign, 
                     pft2biom = PFT2biom, 
                     transf=1)
    colnames(pft$biome) <- "BIOMPOL"
    PFT.output <- cbind(BIOMPOL = data.frame(pft$biome), data.frame(pft$scores))
    Biom.output <- cbind(BIOMPOL = data.frame(pft$biome), data.frame(pft$biomsco))
  }
  
  #### Export results ####
  if(is.null(Save.path) == F){
    saveRDS(PFT.output, gsub("\\.Rds", "_PFT.Rds", Save.path))
    saveRDS(Biom.output, gsub("\\.Rds", "_biomization.Rds", Save.path))}
  
  return(list(PFT.output = PFT.output, Biom.output = Biom.output))
}

Fossil.MAT.prep <- function(MP_fossil, MAge, Type_MAT, Displot, Show.message = F, Save.missing = NULL,
                            Verbose = T, H = 400, W = 700, Save.plot = NULL, Seuil.p = 0.01, Corresp_name){
  #### Settings ####
  if(missing(Displot)){Displot = F}
  if(missing(Corresp_name)){warning("Missing the matrix of the pollen type / transfert function type correspondance.")}
  
  #### Donnees en %TP ou en [C] ####
  MP_fossil[is.na(MP_fossil)] <- 0
  if(unique(apply(MP_fossil, 2, typeof)) == "integer"){
    if(Show.message == T){print("Pollen data in counts.")}
    MP_fossil <- (MP_fossil/rowSums(MP_fossil))
    MP_fossil[MP_fossil <= 0.01] <- 0
    MP_fossil <- (MP_fossil/rowSums(MP_fossil))
  }
  else{
    if(rowSums(MP_fossil[1,]) < 2){
      if(Show.message == T){print("Pollen data in %.")}
      MP_fossil[MP_fossil <= 0.01] <- 0
      MP_fossil <- (MP_fossil/rowSums(MP_fossil))}                    # Fossil data in %TP
    else{              
      if(Show.message == T){print("Pollen data in concentration.")}
      MP_fossil[MP_fossil <= 0.5] <- 0
      MP_fossil <- (MP_fossil/rowSums(MP_fossil))}                    # Fossile pollen [C]
  }
  
  #### Nom taxon correspondant MAT model / fossils samples ####
  a <- Corresp_name[Type_MAT]
  Dup <- Corresp_name[duplicated(Corresp_name$Nom) | duplicated(Corresp_name$Nom, fromLast = T),]
  if(nrow(Dup) > 0){
    Dup <- Dup[order(Dup$Nom),]
    print("**** Be carefull, the following taxa are duplicated within the correspondance table. Please clean them. ****")
    print(Dup)
  }
  
  # print(setdiff(names(MP_fossil),Corresp_name$Nom))
  Keep.taxa <- subset(Corresp_name, a != "ABS")                                                                   # Valeurs Absente de la DB_surface
  names(MP_fossil) <- gsub("\\."," ",names(MP_fossil))                                                            # remplace . en espace
  MP_fossil.cor <- MP_fossil[, which(names(MP_fossil) %in% Keep.taxa$Nom)]                                        # on enleve les taxons absents
  Keep.taxa <- Keep.taxa[match(intersect(names(MP_fossil.cor),Keep.taxa$Nom), Keep.taxa$Nom),]                    # on remplace les noms par ceux de la DB_surf
  if(nrow(Keep.taxa[Keep.taxa[[Type_MAT]] == "",])>0){
    print("**** The following taxa are empty in the FT Type colomn of the correspondance table ! Please fullfil them. ****")
    print(Keep.taxa[Keep.taxa[[Type_MAT]] == "",])
  }
  
  colnames(MP_fossil.cor)<- Keep.taxa[[Type_MAT]]                                                            
  MP_fossil.cor <- as.data.frame(do.call(cbind, by(t(MP_fossil.cor),INDICES=names(MP_fossil.cor),FUN=colSums)))   # on somme les taxons qui appartiennent au meme type
  MP_fossil.cor <- (MP_fossil.cor/rowSums(MP_fossil.cor))
  
  if(length(setdiff(names(MP_fossil), Corresp_name$Nom)) >= 1){
    print(paste("Les taxons suivants ont été enlevés de la base, faute de correspondance dans le fichier :", paste(setdiff(names(MP_fossil), Corresp_name$Nom), collapse = ", " )))
    if(is.null(Save.missing) == F){write.csv(setdiff(names(MP_fossil), Corresp_name$Nom), file = Save.missing)}
    if(Show.message == T){
      print(paste("Les taxons suivants ont été enlevés de la base, car 'ABS' dans le type de FT:", paste(setdiff(names(MP_fossil), Keep.taxa$Nom), collapse = ", " )))
    }
  }
  
  #### Vérif si les échantillons sont bien dans le bon sens (depth ou age) ####
  Common.samples <- intersect(row.names(MP_fossil.cor), row.names(MAge))
  MP_fossil.cor <- MP_fossil.cor[na.omit(match(row.names(MAge), Common.samples)),]
  MAge <- MAge[na.omit(match(row.names(MAge), Common.samples)),]
  
  #### Plot verif ####
  if(Displot == T){
    #### Save plots ####
    if(is.null(Save.plot) == F){
      Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
      dir.create(file.path(Path.to.create), showWarnings = FALSE)
      if(is.null(W) == F & is.null(H) == F){
        pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
      else{pdf(file = Save.plot)}}
    
    #### Plot ####
    mx <- apply(MP_fossil.cor, 2, max)
    MP_fossil.sub <- MP_fossil.cor[, mx > Seuil.p]   # seuil %TP > 5%
    
    Plot.x <- "Age"
    if(any(is.na(MAge[[Plot.x]]))){Plot.x = "Top"}
    rioja::strat.plot(MP_fossil.sub*100, yvar = MAge[[Plot.x]], scale.percent=T, y.rev=T, plot.poly=T, col.poly.line=NA, exag=T, col.exag="auto", col.poly="darkgreen")
    if(is.null(Save.plot) == F){dev.off()}
    
  }
  
  #### Return algo ####
  return(MP_fossil.cor)
}

Fossil.corresp.surface <- function(MP_fossil, MP_surf){
  inter <- intersect(names(MP_fossil), names(MP_surf))
  MP_fossil <- MP_fossil[,inter]
  
  Taxa.miss <- setdiff(names(MP_surf),names(MP_fossil))
  A <- setNames(data.frame(matrix(ncol = length(Taxa.miss), nrow = nrow(MP_fossil))), Taxa.miss) 
  row.names(A) <- row.names(MP_fossil)
  A[is.na(A)] <- 0
  A <- cbind(MP_fossil, A)
  MP_fossil <- A[,sort(names(A))]
  
  return(MP_fossil)
}

FT.core <- function(MCore, MAge, Model.WAPLS = NULL, Model.MAT = NULL, Model.RF = NULL, Model.BRT = NULL,
                    Only.fit = F, LakeName, Select.clim, Fit.val, Nodeparse = F,
                    Ecartype.curve, Model.param.show, Displot = T, Verbose = T, GDGT = F, GDGT.model = NULL,
                    Zone.Clim.span = NULL, Zone.Temp, Save.path, Save.tab, Save.plot = NULL, Save.RDS = F, H = NULL, W = NULL){
  #### Init param ####
  if(missing(Save.tab)){Save.tab = T}
  if(missing(Zone.Temp)){Zone.Temp = rep("U", length(Zone.Clim.span)/2)}
  if(missing(Select.clim)){Select.clim = NULL}
  if(missing(Save.path)){Save.tab = F; Save.path = NULL}
  if(missing(Zone.Clim.span)){Zone.OK = F}
  if(missing(Model.param.show)){Model.param.show = F}
  if(missing(LakeName)){LakeName = "Lake"}
  if(missing(Fit.val)){Fit.val = 0}
  if(missing(Ecartype.curve)){Ecartype.curve = c(F,F,F,F)}
  if(missing(MAge)){MAge = paste(LakeName, seq(1:nrow(MCore)), sep = "_")}
  
  #### Select le model type ####
  if(Nodeparse == F){
    Keep.WAPLS <- deparse(substitute(Model.WAPLS))
    Keep.MAT <- deparse(substitute(Model.MAT))
    Keep.RF <- deparse(substitute(Model.RF))
    Keep.BRT <- deparse(substitute(Model.BRT))
  }
  else{
    if(length(Model.WAPLS) == 0){Model.WAPLS <- NULL}
    else{Keep.WAPLS <- names(Model.WAPLS); Model.WAPLS <- Model.WAPLS[[1]]}
    
    if(length(Model.MAT) == 0){Model.MAT <- NULL; Keep.MAT = NULL}
    else{Keep.MAT <- names(Model.MAT); Model.MAT <- Model.MAT[[1]]}
    
    if(length(Model.BRT) == 0){Model.BRT <- NULL; Keep.BRT = NULL}
    else{Keep.BRT <- names(Model.BRT); Model.BRT <- Model.BRT[[1]]}
    
    if(length(Model.RF) == 0){Model.RF <- NULL; Keep.RF = NULL}
    else{Keep.RF <- names(Model.RF); Model.RF <- Model.RF[[1]]}
  }
  
  if(is.null(Model.BRT) == F){if(any(names(Model.BRT) %in% "Settings") == T){Model.BRT <- Model.BRT[!names(Model.BRT) %in% "Settings"]}}
  if(is.null(Model.RF) == F){if(any(names(Model.RF) %in% c("Settings", "External.CV.Param")) == T){Model.RF <- Model.RF[!names(Model.RF) %in% c("Settings", "External.CV.Param")]}}
  
  if(is.null(Model.WAPLS)== F){Model.type <- Model.WAPLS}
  if(is.null(Model.MAT)== F){Model.type <- Model.MAT}
  if(is.null(Model.RF)== F){Model.type <- Model.RF}
  if(is.null(Model.BRT)== F){Model.type <- Model.BRT}
  
  #### Save plots ####
  if(is.null(Save.plot) == F & Displot == T){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Select param clim ####
  if(is.null(Select.clim) == F){
    Select.clim <- c(Select.clim, "Best.Param")
    Model.type <- Model.type[names(Model.type) %in% Select.clim]
    Model.type$Best.Param <- Model.type$Best.Param[row.names(Model.type$Best.Param) %in% Select.clim,]
    
    if(is.null(Model.BRT) == F){
      Model.BRT <- Model.BRT[names(Model.BRT) %in% Select.clim]
      print(Model.BRT)
      Model.BRT$Best.Param <- Model.BRT$Best.Param[row.names(Model.BRT$Best.Param) %in% Select.clim,]
    }
    
    if(is.null(Model.MAT) == F){
      Model.MAT <- Model.MAT[names(Model.MAT) %in% Select.clim]
      Model.MAT$Best.Param <- Model.MAT$Best.Param[row.names(Model.MAT$Best.Param) %in% Select.clim,]
    }
    
    if(is.null(Model.WAPLS) == F){
      Model.WAPLS <- Model.WAPLS[names(Model.WAPLS) %in% Select.clim]
      Model.WAPLS$Best.Param <- Model.WAPLS$Best.Param[row.names(Model.WAPLS$Best.Param) %in% Select.clim,]
    }
    
    if(is.null(Model.RF) == F){
      Model.RF <- Model.RF[names(Model.RF) %in% Select.clim]
      Model.RF$Best.Param <- Model.RF$Best.Param[row.names(Model.RF$Best.Param) %in% Select.clim,]
    }
  }
  
  #### Graphical settings ####
  Zone.Temp <- rep(Zone.Temp, each=2)
  if(is.null(Zone.Clim.span) == T | length(Zone.Clim.span) != length(Zone.Temp)){
    if(Displot == T){print("Something is wrong with the climate zones. Please check.")}
    Zone.OK = F}
  else{Zone.OK = T}
  
  Tailleplot <-length(Model.type)-1
  if (Tailleplot <= 3){par(mfrow = c(1,Tailleplot))}
  if (Tailleplot == 4){par(mfrow = c(2,2))}
  if (Tailleplot >= 5 & Tailleplot <= 6){par(mfrow = c(2,3))}
  if (Tailleplot >= 7 & Tailleplot <= 9){par(mfrow = c(3,3))}
  Mmodel <- data.frame(Age = MAge)
  MModel.MAT <- data.frame(Age = MAge)
  MModel.RF <- data.frame(Age = MAge)
  MModel.BRT <- data.frame(Age = MAge)
  M.errors.WAPLS <- data.frame(Age = MAge)
  M.errors.MAT <- data.frame(Age = MAge) 
  M.errors.RF <- data.frame(Age = MAge)
  M.errors.BRT <- data.frame(Age = MAge) 
  Full.MAT.RDS <- list()
  Full.WAPLS.RDS <- list()
  Full.RF.RDS <- list()
  Full.BRT.RDS <- list()
  
  #### Loop on the climat param ####
  if(Verbose == F){
    library(lubridate)
    pb = txtProgressBar(min = 1, max = (length(Model.type)-1), width = 40, initial = 0,  style = 3) 
    init <- numeric((length(Model.type)-1))
    end <- numeric((length(Model.type)-1))
  }
  
  print(paste("Prediction for ", LakeName, " with the following models :", Keep.WAPLS, ", ", Keep.MAT, ", ", Keep.BRT, ", ",  Keep.RF, ".", sep = ""))
  for (i in 1:(length(Model.type)-1)){
    if(Verbose == F){init[i] <- Sys.time()}
    LabParamlim <- as.character(gsub("\\p{P}","", deparse(names(Model.type)[i]), perl = TRUE ))
    #### WAPLS ####
    if(is.null(Model.WAPLS) == F){
      NComp.WAPLS = Model.WAPLS$Best.Param[[3]][i]
      if(NComp.WAPLS < 2){NComp.WAPLS <- 2}
      if(Verbose == T){print(paste(round(i/(length(Model.WAPLS)), digits = 2)*100, "% done. The ", LabParamlim, " is modelling with the WAPLS and ", NComp.WAPLS, " parameters.", sep = ""))}
      
      if(mean(rowSums(Model.WAPLS[[i]]$y), na.rm = T) > 1.5){
        MCore.WAPLS <- sqrt(MCore)
        print("Sqrt() transformation of the past data.")}
      else{MCore.WAPLS <- MCore}
      
      Cor.WAPLS = predict(Model.WAPLS[[i]], MCore.WAPLS, npls = NComp.WAPLS, sse = T, nboot = 1000, verbose = F)
      Mmodel[i+1] <- cbind(Cor.WAPLS$fit[,NComp.WAPLS])
      colnames(Mmodel)[i+1] <- LabParamlim
      M.errors.WAPLS[i+1] <- cbind(Cor.WAPLS$SEP.boot[,NComp.WAPLS])
      colnames(M.errors.WAPLS)[i+1] <- LabParamlim
      Full.WAPLS.RDS[[i]] <- Cor.WAPLS
      names(Full.WAPLS.RDS)[[i]] <- LabParamlim}
    
    #### MAT ####
    if(is.null(Model.MAT) == F){
      NComp.MAT = Model.MAT$Best.Param[[3]][i]
      if(NComp.MAT < 4){NComp.MAT <- 4}
      if(Verbose == T){print(paste("The ", LabParamlim, " is predicted with the MAT and ", NComp.MAT, " analogues.", sep = ""))}
      
      if(mean(rowSums(Model.MAT[[i]]$y), na.rm = T) > 1.5){
        MCore.MAP <- sqrt(MCore)
        print("Sqrt() transformation of the past data.")}
      else{MCore.MAP <- MCore}
      
      Cor.MAT = predict(Model.MAT[[i]], MCore.MAP, k = NComp.MAT, sse = T, nboot = 1000, verbose = F)
      MModel.MAT[i+1] <- cbind(Cor.MAT$fit[,2])
      M.errors.MAT[i+1] <- cbind(Cor.MAT$SEP.boot[,2])
      colnames(MModel.MAT)[i+1] <- LabParamlim
      colnames(M.errors.MAT)[i+1] <- LabParamlim
      Full.MAT.RDS[[i]] <- Cor.MAT
      names(Full.MAT.RDS)[[i]] <- LabParamlim}
    
    #### Random forest (RF) ####
    if(is.null(Model.RF) == F){
      if(Verbose == T){print(paste("The ", LabParamlim, " is predicted with the RF.", sep = ""))}
      Cor.RF = predict(Model.RF[[i]], MCore, verbose = F)
      MSE <- Model.RF[[i]]$mse
      RMSE.RF <- sqrt(MSE[length(MSE)])
      MModel.RF <- cbind(MModel.RF, Cor.RF)
      colnames(MModel.RF)[i+1] <- LabParamlim
      Full.RF.RDS[[i]] <- Cor.RF
      names(Full.RF.RDS)[[i]] <- LabParamlim}
    
    #### Boosted Regression Tree (BRT) ####
    if(is.null(Model.BRT) == F){
      if(Verbose == T){print(paste("The ", LabParamlim, " is predicted with the BRT.", sep = ""))}
      MCore.i <- MCore[, names(MCore) %in% colnames(Model.BRT[[i]]$data$x.order)]
      Cor.BRT <- gbm::predict.gbm(Model.BRT[[i]], MCore.i, n.trees = Model.BRT[[i]]$gbm.call$best.trees, type="response", sse = T, verbose = F)
      RMSE.BRT <- Model.BRT$Best.Param[i,2]
      MModel.BRT<- cbind(MModel.BRT, Cor.BRT)
      colnames(MModel.BRT)[i+1] <- LabParamlim
      Full.BRT.RDS[[i]] <- Cor.BRT
      names(Full.BRT.RDS)[[i]] <- LabParamlim
    }
    
    #### Plot résultats graphiques ####
    if(Displot == T){
      #### Val Min / Max ####
      ymin = 0
      ymax = 0
      if(Ecartype.curve[1] == F){
        
        if(is.null(Model.WAPLS) == F){
          ymin = min(Cor.WAPLS$fit[,NComp.WAPLS], na.rm = T)
          ymax = max(Cor.WAPLS$fit[,NComp.WAPLS], na.rm = T)}
        
        if(is.null(Model.MAT) == F & Ecartype.curve[2] == F){
          ymin = min(ymin, Cor.MAT$fit[,2], na.rm = T)
          ymax = max(ymax, Cor.MAT$fit[,2], na.rm = T)}
        
        if(is.null(Model.RF) == F & Ecartype.curve[3] == F){
          ymin = min(ymin, min(Cor.RF, na.rm = T), na.rm = T)
          ymax = max(ymax, max(Cor.RF, na.rm = T), na.rm = T)}
        
        if(is.null(Model.BRT) == F & Ecartype.curve[4] == F){
          ymin = min(ymin, min(Cor.BRT, na.rm = T), na.rm = T)
          ymax = max(ymax, max(Cor.BRT, na.rm = T), na.rm = T)}
        
        fullY <- abs(ymax) - abs(ymin)
        ymax = ymax + 0.05*fullY
      }
      else{
        if(is.null(Model.WAPLS) == F){
          ymin = min(Cor.WAPLS$fit[,NComp.WAPLS] - Cor.WAPLS$SEP.boot[,1], na.rm = T)
          ymax = max(Cor.WAPLS$fit[,NComp.WAPLS] + Cor.WAPLS$SEP.boot[,1], na.rm = T)}
        
        if(is.null(Model.MAT) == F){
          ymin = min(ymin, na.omit(Cor.MAT$fit[,2] - Cor.MAT$SEP.boot[,1]), na.rm = T)
          ymax = max(ymax, na.omit(Cor.MAT$fit[,2] + Cor.MAT$SEP.boot[,1]), na.rm = T)}
        
        if(is.null(Model.RF) == F){
          ymin = min(ymin, min(Cor.RF - RMSE.RF, na.rm = T), na.rm = T)
          ymax = max(ymax, max(Cor.RF + RMSE.RF, na.rm = T), na.rm = T)}
        
        if(is.null(Model.BRT) == F){
          ymin = min(ymin, c(Cor.BRT - RMSE.BRT), na.rm = T)
          ymax = max(ymax, max(Cor.BRT + RMSE.BRT, na.rm = T), na.rm = T)}
        
        fullY <- abs(ymax) - abs(ymin)
        ymax = ymax + 0.05*fullY
      }
      
      #### Plot Mean value ####
      if(is.null(Model.WAPLS) == F & Only.fit == F){Y <- Cor.WAPLS$fit[,NComp.WAPLS]}
      else{Y <- rep(NA, length(MAge))}
      
      plot(MAge, Y, 
           xlim = c(min(MAge),max(MAge)), 
           ylim = c(ymin, ymax), 
           type = "l", 
           ylab = LabParamlim, 
           xlab = "Time (yr cal BP)", 
           col = 1, 
           las = 0, lwd=1, bty="n")
      
      #### Legend ####
      if(Model.param.show == T){
        fullT <- abs(min(MAge)) + abs(max(MAge))
        x1 = 0.3*fullT
        x2 = 0.4*fullT
        x4 = 0.55*fullT
        x5 = 0.65*fullT
        
        ymax.lab1 <- ymax 
        ymax.lab2 <- ymax - 0.05*ymax
        
        #### Legend WAPLS ####
        if(is.null(Model.WAPLS) == F){
          mylabel1 = Keep.WAPLS
          mylabel2 = bquote(npls == .(Model.WAPLS[[1]][["npls"]]) ~ "," ~
                              k == .(Model.WAPLS[[length(Model.WAPLS)]][i,3]) ~ "," ~
                              italic(R)^2 == .(format(Model.WAPLS[[length(Model.WAPLS)]][i,4], digits = 2)) ~ "," ~
                              RMSE == .(format(Model.WAPLS[[length(Model.WAPLS)]][i,5], digits = 2)))
          
          text(x = x1, y = ymax.lab1, labels = mylabel1, font = 2, cex = 0.8, pos = 1)
          text(x = x2, y = ymax.lab1, labels = mylabel2, cex = 0.8, pos = 1)
        }
        
        #### Legend MAT ####
        if(is.null(Model.MAT) == F){
          mylabel1b = Keep.MAT
          mylabel3b = bquote(k == .(Model.MAT[[length(Model.MAT)]][i,3]) ~ "," ~
                               italic(R)^2 == .(format(Model.MAT[[length(Model.MAT)]][i,4], digits = 2)) ~ "," ~
                               RMSE == .(format(Model.MAT[[length(Model.MAT)]][i,5], digits = 2)))
          
          
          text(x = x1, y = ymax.lab2, labels = mylabel1b, col = "royalblue", font = 2, cex = 0.8)
          text(x = x2, y = ymax.lab2, labels = mylabel3b, cex = 0.8)
        }
        
        #### Legend RF ####
        if(is.null(Model.RF) == F){
          mylabel1b = Keep.RF
          mylabel2b = bquote(
            italic(R)^2 == .(format(Model.RF[[length(Model.RF)]][i,1], digits = 2)) ~ "," ~
              RMSE == .(format(Model.RF[[length(Model.RF)]][i,2], digits = 2)))
          
          
          text(x = x4, y = ymax.lab2, labels = mylabel1b, col = "darkorange", font = 2, cex = 0.8, pos = 1)
          text(x = x5, y = ymax.lab2, labels = mylabel2b, cex = 0.8, pos = 1)
        }
        
        
        #### Legend BRT ####
        if(is.null(Model.BRT) == F){
          mylabel1b = Keep.BRT
          mylabel2b = bquote(Nb.tree == .(Model.BRT[[i]][["call"]][["ntree"]]) ~ "," ~
                               italic(R)^2 == .(format(Model.BRT[[length(Model.BRT)]][i,1], digits = 2)) ~ "," ~
                               RMSE == .(format(Model.BRT[[length(Model.BRT)]][i,2], digits = 2)))
          
          
          text(x = x4, y = ymax.lab1, labels = mylabel1b, col = "darkgreen", font = 2, cex = 0.8, pos = 1)
          text(x = x5, y = ymax.lab1, labels = mylabel2b, cex = 0.8, pos = 1)
        }
      }
      
      #### Plot Model MAT add ####
      if (is.null(Model.MAT) == F & Only.fit == F){
        lines(MAge, Cor.MAT$fit[,2], col = "royalblue", lwd=1, las=0)}
      
      #### Plot Model RF add ####
      if(is.null(Model.RF) == F & Only.fit == F){
        lines(MModel.RF[[1]], Cor.RF, col = "darkorange", lwd=1, las=0)}
      
      #### Plot Model BRT add ####
      if(is.null(Model.BRT) == F & Only.fit == F){
        lines(MModel.BRT[[1]], Cor.BRT, col = "darkgreen", lwd=1, las=0)}
      
      #### Plot fitting ####
      if(Fit.val > 0){
        if (is.null(Model.WAPLS) == F){
          Curve.fit = lowess(Cor.WAPLS$fit[,NComp.WAPLS], f = Fit.val)
          lines(MAge, Curve.fit$y, col=1, lwd=2)}
        if (is.null(Model.MAT) == F){
          Curve.fit.MAT = lowess(Cor.MAT$fit[,2], f = Fit.val)
          lines(MAge, Curve.fit.MAT$y, col= "royalblue", lwd=2)}
        if (is.null(Model.RF) == F){
          Curve.fit.MAT = lowess(Cor.RF, f = Fit.val)
          lines(MAge, Curve.fit.MAT$y, col = "darkorange", lwd=2)}
        if (is.null(Model.BRT) == F){
          Curve.fit.MAT = lowess(Cor.BRT, f = Fit.val)
          lines(MAge, Curve.fit.MAT$y, col = "darkgreen", lwd=2)}
      }
      
      #### Interval WAPLS ####
      if (Ecartype.curve[1] == T & is.null(Model.WAPLS) == F){
        lines(MAge, Cor.WAPLS$fit[,NComp.WAPLS] + Cor.WAPLS$SEP.boot[,1], lwd=.4, lty = "dashed")
        lines(MAge, Cor.WAPLS$fit[,NComp.WAPLS] - Cor.WAPLS$SEP.boot[,1], lwd=.4, lty = "dashed")
      }
      
      #### Interval MAT ####
      if (Ecartype.curve[2] == T & is.null(Model.MAT) == F){
        lines(MAge, Cor.MAT$fit[,2] + Cor.MAT$SEP.boot[,1], lwd=.4, col = "royalblue", lty = "dashed")
        lines(MAge, Cor.MAT$fit[,2] - Cor.MAT$SEP.boot[,1], lwd=.4, col = "royalblue", lty = "dashed")
      }
      
      #### Interval RF ####
      if (Ecartype.curve[3] == T & is.null(Model.RF) == F){
        lines(MAge, Cor.RF + RMSE.RF, lwd=.4, col = "darkorange", lty = "dashed")
        lines(MAge, Cor.RF - RMSE.RF, lwd=.4, col = "darkorange", lty = "dashed")
      }
      
      #### Interval BRT ####
      if (Ecartype.curve[4] == T & is.null(Model.BRT) == F){
        lines(MAge, Cor.BRT + RMSE.BRT, lwd=.4, col = "darkgreen", lty = "dashed")
        lines(MAge, Cor.BRT - RMSE.BRT, lwd=.4, col = "darkgreen", lty = "dashed")
      }
      
      #### Plot climate zones ####
      if (Zone.OK == T){
        for(j in 1:length(Zone.Clim.span)){
          if(as.logical(j%%2) == T){        # seulement les j impairs
            a <- which(MAge>Zone.Clim.span[j])
            b <- which(MAge>Zone.Clim.span[j+1])
            chaud = rgb(1, 0, 0, 0.08)
            froid = rgb(0, 0, 1, 0.08)
            unknow = rgb(0.3, 0.3, 0.3, 0.08)
            if(Zone.Temp[j]=="C"){colTemp <- froid}
            if(Zone.Temp[j]=="W"){colTemp <- chaud}
            if(Zone.Temp[j]=="U"){colTemp <- unknow}
            polygon(c(Zone.Clim.span[j],Zone.Clim.span[j+1], Zone.Clim.span[j+1],Zone.Clim.span[j]), c(ymin, ymin, ymax, ymax), border = NA, col = colTemp)
          }}}
    }
    if(Verbose == F){
      end[i] <- Sys.time()
      setTxtProgressBar(pb, i)
      time <- round(seconds_to_period(sum(end - init)), 0)
      est <- (length(Model.type)-1) * (mean(end[end != 0] - init[init != 0])) - time
      remainining <- round(seconds_to_period(est), 0)
      cat(paste(" // Execution time:", time,
                " // Estimated time remaining:", remainining), "")}
    
  }
  
  if(Verbose == F){close(pb);library(beepr)} 
  
  #### Save / export data ####
  row.names(Mmodel) <- row.names(MCore)
  par(mfrow = c(1,1))
  if(Save.tab == T){
    Path.to.create.csv <- gsub("(.*/).*\\.csv.*","\\1", Save.path)
    dir.create(file.path(Path.to.create.csv), showWarnings = FALSE)
    
    if(is.null(Model.WAPLS) == F){
      #### Save WAPLS ####
      Save.Model.WAPLS <- gsub("\\.", "_", Keep.WAPLS)
      add.to.path <- paste("_", Save.Model.WAPLS, ".csv", sep = "")
      add.to.path.error <- paste("_SEP_", Save.Model.WAPLS, ".csv", sep = "")
      Save.path1 <- gsub("\\.csv", add.to.path, Save.path)
      Save.path.error <- gsub("\\.csv", add.to.path.error, Save.path)
      write.table(Mmodel, file = Save.path1, row.names=T, col.names=NA, sep=",", dec = ".")
      write.table(M.errors.WAPLS, file = Save.path.error, row.names=T, col.names=NA, sep=",", dec = ".")
    }}
  
  if(is.null(Model.MAT) == F){
    #### Save MAT ####
    if(Save.tab == T){
      Save.Model.WAPLS <- gsub("\\.", "_", Keep.MAT)
      add.to.path <- paste("_", Save.Model.WAPLS, ".csv", sep = "")
      add.to.path.error <- paste("_SEP_", Save.Model.WAPLS, ".csv", sep = "")
      Save.path2 <- gsub("\\.csv", add.to.path, Save.path)
      Save.path2.error <- gsub("\\.csv", add.to.path.error, Save.path)
      write.table(MModel.MAT, file = Save.path2, row.names=T, col.names=NA, sep=",", dec = ".")
      write.table(M.errors.MAT, file = Save.path2.error, row.names=T, col.names=NA, sep=",", dec = ".")
    }
    
    #### Add MAT to function return ####
    Mmodel = list(Mmodel, MModel.MAT, M.errors.WAPLS, M.errors.MAT)
    Save.DB.name <- gsub(".*\\.","", Keep.WAPLS)
    labtot <- c(Keep.WAPLS, Keep.MAT, paste("SEP.WAPLS", Save.DB.name, sep = "."), paste("SEP.MAT", Save.DB.name, sep = "."))
    names(Mmodel) <- labtot}
  
  if(is.null(Model.RF) == F){
    #### Save RF ####
    if(Save.tab == T){
      Save.Model.WAPLS <- gsub("\\.", "_", Keep.RF)
      add.to.path <- paste("_", Save.Model.WAPLS, ".csv", sep = "")
      Save.path2 <- gsub("\\.csv", add.to.path, Save.path)
      write.table(MModel.RF, file = Save.path2, row.names=T, col.names=NA, sep=",", dec = ".")
    }
    
    #### Add RF to function return ####
    Mmodel <- append(Mmodel, list(MModel.RF))
    names(Mmodel)[length(Mmodel)] <- Keep.RF
  }
  
  if(is.null(Model.BRT) == F){
    #### Save BRT ####
    if(Save.tab == T){
      Save.Model.BRT <- gsub("\\.", "_", Keep.BRT)
      add.to.path <- paste("_", Save.Model.BRT, ".csv", sep = "")
      Save.path2 <- gsub("\\.csv", add.to.path, Save.path)
      write.table(MModel.BRT, file = Save.path2, row.names=T, col.names = NA, sep=",", dec = ".")
    }
    
    #### Add BRT to function return ####
    Mmodel <- append(Mmodel, list(MModel.BRT))
    names(Mmodel)[length(Mmodel)] <- Keep.BRT
  }
  
  Total.model <- list(WAPLS = Full.WAPLS.RDS, MAT = Full.MAT.RDS, BRT = Full.BRT.RDS, RF = Full.RF.RDS)
  
  #### Save format RDS ####
  if(Save.RDS == T & is.null(Save.path) == F){
    if(GDGT == F & exists("Save.DB.name") == F){Save.DB.name <- gsub(".*\\.","", Keep.WAPLS)}
    if(GDGT == T){
      if(is.null(Model.BRT) == F){Save.DB.name <- gsub(".*\\.","", Keep.BRT)}
      if(is.null(Model.WAPLS) == F){Save.DB.name <- gsub(".*\\.","", Keep.WAPLS)}
      if(is.null(Model.MAT) == F){Save.DB.name <- gsub(".*\\.","", Keep.MAT)}
      if(is.null(Model.RF) == F){Save.DB.name <- gsub(".*\\.","", Keep.RF)}
      
      if(is.null(GDGT.model) == F){Save.DB.name <- paste(GDGT.model, Save.DB.name, sep = "_")} 
    }
    
    Save.path.RDS <- paste(gsub("\\.csv", paste("_", Save.DB.name, sep = ""), Save.path), ".Rds", sep = "")
    Save.path.RDS.full <- paste(gsub("\\.csv", paste("_", Save.DB.name, "_full", sep = ""), Save.path), ".Rds", sep = "")
    Path.to.create2 <- gsub("(.*/).*\\.Rds.*","\\1", Save.path.RDS)
    dir.create(file.path(Path.to.create2), showWarnings = F)
    
    saveRDS(Total.model, Save.path.RDS.full)
    saveRDS(Mmodel, Save.path.RDS)}
  
  #### End ####
  if(is.null(Save.plot) == F){dev.off()}
  if(GDGT == F){return(Total.model)}
  else{return(Mmodel)}
}

Cal.Mtern <- function(M){
  M <- M[,grep("^f.I", colnames(M))]
  III <- grep("III", names(M))
  III.II <- grep("II", names(M))
  I.II.III <- grep("I", names(M))
  II <- setdiff(III.II, III)  
  I <- setdiff(I.II.III, III.II)
  Mtern <- M[0]
  Mtern[["P.hexa"]] <- rowSums(M[,III])
  Mtern[["P.penta"]] <- rowSums(M[,II])
  Mtern[["P.tetra"]] <- rowSums(M[,I])
  Mtern <- data.frame(apply(Mtern, 2, function(x) x/rowSums(Mtern))) 
  return(Mtern)}

Add.ML.to.cores <- function(Mcore, Mmin = NULL, Mmax = NULL, Mcore.ML, Keep.param = NULL, Msurf.ML = NULL) {
  if(is.null(Keep.param) == T){
    Remove <- c("Age", "Depth", "Top", "Bottom")
    Keep.param <- names(Mcore.ML[[1]])
    Keep.param <- Keep.param[!Keep.param %in% Remove]
  }
  
  for(i in 1:length(Mcore.ML)) {
    Mcore.ML.i <- Mcore.ML[[i]]
    if(is.null(Msurf.ML) == F){Msurf.ML.i <- Msurf.ML[[i]]}
    
    Lab.i <- names(Mcore.ML)[i]
    Mcore.ML.i <- Mcore.ML.i[c(which(names(Mcore.ML.i) %in% Keep.param))]
    names(Mcore.ML.i) <- paste(names(Mcore.ML.i), Lab.i, sep = "_")
    
    Mcore <- cbind(Mcore, Mcore.ML.i)
    if(is.null(Mmin) == F){
      if(identical(names(Mcore.ML), names(Msurf.ML)) == F){print("Attention, les matrices Msurf et Mcore doivent être de même taille et porter les mêmes noms!")}
      else{
        Mcore.ML.i.MaxI <- Mcore.ML.i
        Mcore.ML.i.MinI <- Mcore.ML.i
        for(j in 1:length(Keep.param)){
          Mcore.ML.i.MaxI[j] <- Mcore.ML.i.MaxI[j] + Msurf.ML.i$Best.Param$RMSE.choice[which(Keep.param[j] == row.names(Msurf.ML.i$Best.Param))]/2
          Mcore.ML.i.MinI[j] <- Mcore.ML.i.MinI[j] - Msurf.ML.i$Best.Param$RMSE.choice[which(Keep.param[j] == row.names(Msurf.ML.i$Best.Param))]/2
        }
        Mmin <- cbind(Mmin, Mcore.ML.i.MaxI)
        Mmax <- cbind(Mmax, Mcore.ML.i.MinI)
      }
    }
  }
  
  Mexp <- list(Mcore)
  if(is.null(Mmin) == F){Mexp[[length(Mexp) + 1]] <- Mmin}
  if(is.null(Mmax) == F){Mexp[[length(Mexp) + 1]] <- Mmax}
  
  return(Mexp)
}

GDGT.local.rescale <- function(Keep.models = "BRT", List.params = c("MAAT"), Msurf = NULL, Mclim = NULL, Mpaleo = NULL, Mpaleo.MaxI = NULL, Mpaleo.MinI = NULL){
  for(i in 1:length(Keep.models)){
    for(j in 1:length(List.params)){
      To.correct <- names(Msurf)[grepl(List.params[j], names(Msurf))&grepl(Keep.models[i], names(Msurf))]
      Clim.j <- List.params[j]
      for(k in 1:length(To.correct)){
        if(length(To.correct[k]) > 0){
          if(is.na(To.correct[k]) == F){
            Mcor <- cbind(Msurf[To.correct[k]], Mclim[Clim.j])
            reg <- lm(Mcor[[2]] ~ Mcor[[1]], data = Mcor)
            Mcor[3] <- Mcor[To.correct[k]]*reg$coefficients[[2]] + reg$coefficients[[1]]
            Msurf[paste(To.correct[k], "LS", sep = "_")] <- Mcor[3]
            if(is.null(Mpaleo) == F){Mpaleo[paste(To.correct[k], "LS", sep = "_")] <- Mpaleo[To.correct[k]]*reg$coefficients[[2]] + reg$coefficients[[1]]}
            if(is.null(Mpaleo.MaxI) == F){Mpaleo.MaxI[paste(To.correct[k], "LS", sep = "_")] <- Mpaleo.MaxI[To.correct[k]]*reg$coefficients[[2]] + reg$coefficients[[1]]}
            if(is.null(Mpaleo.MinI) == F){Mpaleo.MinI[paste(To.correct[k], "LS", sep = "_")] <- Mpaleo.MinI[To.correct[k]]*reg$coefficients[[2]] + reg$coefficients[[1]]}
          }
        }
      }
    }
  }
  
  Mexport <- list(Msurf)
  if(is.null(Mpaleo) == F){Mexport[[2]] <- Mpaleo}
  if(is.null(Mpaleo.MaxI) == F){Mexport[[3]] <- Mpaleo.MaxI}
  if(is.null(Mpaleo.MinI) == F){Mexport[[4]] <- Mpaleo.MinI}
  return(Mexport)
}

Combine.ML.cluster <- function(Cluster.prediction, List.models, Model.lab, GDGT.paleo, Plot.y = "Age", Param.clim = "MAAT", 
                               Highlight.combined = F, Save.path = NULL, Method = "BRT", return.plot = F, Panel.annot = NULL, Time.res = NULL,
                               Compare.curve = NULL, Cluster.prob = "Both", Time.lim = NULL, Surf.val = NULL, Dot.size = 1.5, Time.in.k = F,
                               Core.name = NULL, Plot.y.lab = Plot.y, Show.proba = T, Facet = T, Only.best = T, H = 900, W = 500, Save.plot = NULL){
  #### Cluster predictions ####
  Br.GDGT.paleo <- GDGT.paleo[grepl("f.I", names(GDGT.paleo)) & !grepl("_7Me", names(GDGT.paleo))]
  Br.GDGT.paleo <- data.frame(t(Br.GDGT.paleo))
  Br.GDGT.paleo <- apply(Br.GDGT.paleo, 2, MESS::round_percent)
  Br.GDGT.paleo <- data.frame(t(Br.GDGT.paleo/100))
  RF_class <- stats::predict(Cluster.prediction, Br.GDGT.paleo)
  RF_class_prob <- stats::predict(Cluster.prediction, Br.GDGT.paleo, type = "prob")
  Br.GDGT.paleo$Pred.cluster <- RF_class
  Br.GDGT.paleo$Plot.y <- GDGT.paleo[[Plot.y]]
  
  #### Matrix full ####
  List.models <- Map(function(df, Model){df$Model <- Model; return(df)}, List.models, Model.lab)
  M <- do.call(rbind, List.models)
  All.param <- setdiff(names(M), c(Plot.y, "Model"))
  M$Pred.cluster <- Br.GDGT.paleo$Pred.cluster[match(M[[Plot.y]], Br.GDGT.paleo$Plot.y)]
  M$Model <- factor(M$Model, ordered = T, levels = unique(M$Model))
  
  #### Matrix full (weighted) ####
  List.models[[1]] <- List.models[[1]][c(Plot.y, Param.clim)]
  Mw <- cbind(List.models[[1]], MAAT.karid = List.models[[2]][[Param.clim]], MAAT.kwet = List.models[[3]][[Param.clim]], RF_class_prob)
  names(Mw)[c(1:2)] <- c("Plot.y", "Param.clim")
  Mw$Combine.weighted <- Mw$MAAT.kwet*Mw$`K-cold/wet` + Mw$MAAT.karid*Mw$`K-warm/arid`
  Mw$Pred.cluster <- ifelse(Mw$`K-warm/arid` >= 0.5, "K-warm/arid", "K-cold/wet")
  Mw$Model <- "Combine-Weighted"
  Mw$Param.clim <- Mw$Combine.weighted
  Keep <- c(Plot.y, Param.clim, "Model", "Pred.cluster")
  M <- M[Keep]
  names(M)[c(1:2)] <- c("Plot.y", "Param.clim")
  M <- dplyr::full_join(M, Mw[c(1,2,9,8)], by = join_by(Plot.y, Param.clim, Model, Pred.cluster))
  
  if(Only.best == T){
    M <- M[M$Model  %in% c("ACADB", "Combine-Weighted"),]  
  }
  M$Model <- factor(M$Model, ordered = T, levels = unique(M$Model))
  
  RF_class_prob <- data.frame(RF_class_prob)
  RF_class_prob$Plot.y <- GDGT.paleo[[Plot.y]]
  
  #### Graphical settings ####
  My_colors <- c("Combined" = "grey20", "Combine-Weighted" = "darkred", 
                 "ACADB" = "bisque3", "K-cold/wet" = "royalblue", "K-warm/arid" = "darkorange", "MAAT_Chen_Tjk" = "#6b00b2ff", 
                 "MAAT_soil_Naaf" = "darkgreen", "MAAT_LSun" = "#6b00b2ff", "MAF_MSosa" = "#6b00b2ff", "MAF_meth_Raberg" = "darkgreen", "MAF_full_Raberg" = "aquamarine2",
                 "MAAT_mr_DJ" = "#d1b336ff", "MAAT_DJ_5Me" = "darkolivegreen3", "MAAT_NMSDB_mr5" = "aquamarine2")
  
  My_labs <- c("Combined" = "Combined", "Combine-Weighted" = "BRT(combined)", 
               "MAAT_soil_Naaf" = "MBT'5Me (Naafs et al., 2017)", "MAAT_LSun" = "MBT/CBT lake (Sun et al., 2011)",
               "ACADB" = paste(Method, "(ACADB)"), 
               "K-cold/wet" = paste(Method, "(K-cold/wet)"),
               "K-warm/arid" = paste(Method, "(K-warm/arid)"),
               "MAAT_mr_DJ" = "MR (De Jonge et al., 2014)", "MAAT_DJ_5Me" = "MBT'5Me (De Jonge et al., 2014)", "MAAT_NMSDB_mr5" = "MR mong.")
  
  if(Facet == T){My_facet <- facet_geochem_grid(vars(Model))}
  else{My_facet <- NULL; H <- H/2}
  
  if(is.null(Plot.y.lab) == F){
    if(is.null(Compare.curve) == F){
      Age.lab.2 <- Plot.y.lab; Age.lab.1 <- NULL}
    else{Age.lab.1 <- Plot.y.lab}
    
  }
  else{
    if(is.null(Compare.curve) == F){
      Age.lab.2 <- Plot.y; Age.lab.1 <- NULL}
    else{Age.lab.1 <- Plot.y}
  }
  
  if(Param.clim %in% c("MAAT", "MAF")){
    Clim.lab.1 <- paste(Param.clim, " (°C)\n", Method, sep = "")
    Clim.lab.2 <- paste(Param.clim, " (°C)\ncomparison", sep = "")
  }
  else{
    Clim.lab.1 <- paste(Param.clim, "\n", Method, sep = "")
    Clim.lab.2 <- paste(Param.clim, "\ncomparison", sep = "")
  }
  
  if(is.null(Surf.val) == F){
    Surf.line <- geom_hline(yintercept = Surf.val, linetype = "dotdash", linewidth = 0.55, color = "black", alpha = 0.55)
  }
  else{Surf.line <- NULL}
  
  if(is.null(Core.name) == F){
    My_title <- ggtitle(Core.name)
  }
  else{My_title <- NULL}
  
  if(Highlight.combined == T){
    Double.line <- geom_line(data = Mw, aes(x = Plot.y, y = Combine.weighted, group = Model, color = Model), linewidth = .8)}
  else{Double.line <- NULL}
  
  if(is.null(Time.lim) == F){Xlim <- xlim(Time.lim)}
  else{Xlim <- NULL}
  
  if(is.null(Plot.y.lab) == T){
    X.title <- element_blank()
  }
  else{X.title <- element_text()}
  
  if(is.null(Panel.annot) == F){
    if(length(Panel.annot) == 1){
      Annot.in.1 <- ggplot2::annotate("text", x = Inf, y = Inf, label = paste("(", Panel.annot, "1)", sep = ""), hjust = 1.1, vjust = 1.5, size = 4.5)
      Annot.in.2 <- ggplot2::annotate("text", x = Inf, y = Inf, label = paste("(", Panel.annot, "2)", sep = ""), hjust = 1.1, vjust = 1.5, size = 4.5)
      Annot.in.3 <- ggplot2::annotate("text", x = Inf, y = Inf, label = paste("(", Panel.annot, "3)", sep = ""), hjust = 1.1, vjust = 1.5, size = 4.5)
      
      if(Show.proba == T){Annot.in.cp <- Annot.in.1; Annot.in.main <- Annot.in.2; Annot.in.compar <- Annot.in.3}
      else{Annot.in.main <- Annot.in.1; Annot.in.compar <- Annot.in.2}
    }
    else{print("A faire !!!!")}
  }
  else{Annot.in.cp <- NULL; Annot.in.main <- NULL; Annot.in.compar <- NULL}
  
  if(is.null(Time.res) == F & is.null(Time.lim) == F){
    if(Time.in.k == T){
      My_time_breaks <- paste(seq(Time.lim[1]/1000, Time.lim[2]/1000, Time.res/1000), "k", sep = "")}
    else{My_time_breaks <- seq(Time.lim[1], Time.lim[2], Time.res)}
    My_scale_time <- scale_x_continuous(breaks = seq(Time.lim[1], Time.lim[2], Time.res), labels = My_time_breaks, limits = Time.lim)}
  else{My_scale_time <- NULL}
  
  #### Comparative curves ####
  if(is.null(Compare.curve) == F){
    Mc <- GDGT.paleo[c(Plot.y, Compare.curve)]
    Mc$'Combine-Weighted' <- Mw$Combine.weighted
    Mc <- melt(Mc, id = Plot.y)
    Ticks.1 <- element_blank(); Text.1 <- element_blank()
    
    Mc$Pred.cluster <- Br.GDGT.paleo$Pred.cluster[match(Mc[[Plot.y]], Br.GDGT.paleo$Plot.y)]
    names(Mc)[c(1:2)] <- c("Plot.y", "Model")
    Mc$Model <- factor(Mc$Model, ordered = T, levels = unique(Mc$Model))
    p.comp <- ggplot(Mc, aes(x = Plot.y, y = value, group = Model))+
      My_facet+ Surf.line+ Xlim+ Annot.in.compar + 
      geom_point(aes(color = Model), show.legend = T, size = Dot.size, shape = 16)+
      geom_line(aes(color = Model))+ xlab(Age.lab.2)+ ylab(Clim.lab.2)+
      Double.line +
      My_scale_time +
      scale_color_manual(values = My_colors, label = My_labs, name = NULL)
  }
  else{Ticks.1 <- element_line(); Text.1 <- element_text()}
  
  #### Plot (main) ####
  p.main <- ggplot(M, aes(x = Plot.y, y = Param.clim, group = Model))+
    My_facet+ Surf.line+ Xlim+ Annot.in.main +
    geom_point(aes(color = Pred.cluster), shape = 16, size = Dot.size)+
    geom_line(aes(color = Model))+ 
    Double.line+
    My_scale_time +
    
    xlab(Age.lab.1)+ ylab(Clim.lab.1)+
    scale_color_manual(values = My_colors, label = My_labs, name = NULL)+
    theme(axis.text.x = Text.1, axis.ticks.x = Ticks.1)
  
  #### Plot (proba. cluster) ####
  if(Show.proba == T){
    RF_class_prob2 <- RF_class_prob
    if(Cluster.prob == "K-warm/arid"){
      RF_class_prob2$Prob.to.show <- RF_class_prob2$K.warm.arid
      Estim.col <- "darkorange"; My_col_scale <- NULL
      My_area <- geom_area(fill = Estim.col, position = "stack")}
    if(Cluster.prob == "K-cold/wet"){
      RF_class_prob2$Prob.to.show <- RF_class_prob2$K.cold.wet
      Estim.col <- "royalblue"; My_col_scale <- NULL
      My_area <- geom_area(fill = Estim.col, position = "stack")}
    if(Cluster.prob == "Both"){
      RF_class_prob2 <- melt(RF_class_prob2, id = "Plot.y")
      names(RF_class_prob2)[names(RF_class_prob2) == "value"] <- "Prob.to.show"
      My_area <- geom_area(aes(fill = variable), position = "stack", color = NA, alpha = .7)
      My_col <- c("K.warm.arid" = "darkorange",
                  "K.cold.wet" = "royalblue")
      My_col_scale <- scale_fill_manual(values = My_col, drop = T)
    }
    
    RF_class_prob2$Prob.to.show <- round(RF_class_prob2$Prob.to.show*100, digits = 0)
    p.prob <- ggplot(RF_class_prob2, aes(x = Plot.y, y = Prob.to.show))+
      My_area+ My_title+ Xlim+ Annot.in.cp +
      My_col_scale+
      scale_y_continuous(limits = c(0,100), breaks = c(0,50,100))+
      xlab(NULL)+ylab("Cluster\nprob.")+
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
    
    p.main <- p.prob / p.main
  }
  
  #### Merge plots ####
  if(is.null(Compare.curve) == F){
    p.main <- p.main / p.comp + plot_layout(guides = "collect")
  }
  
  p.main <- p.main&
    theme(panel.background = element_rect(fill = NA, color = "black"), legend.key = element_rect(fill = NA, color = NA),
          panel.grid = element_blank(), axis.title.x = X.title, 
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'pt'),
          legend.position = "bottom")&
    guides(color = guide_legend(nrow = 3))
  
  #### Generalized combined plot ####
  Mexp <- data.frame(List.models[[1]][Plot.y])
  for(i in 1:length(All.param)){
    Clim.i <- All.param[i]
    Mw <- cbind(List.models[[1]], M2 = List.models[[2]][[Clim.i]], M1 = List.models[[3]][[Clim.i]], RF_class_prob)
    
    Mexp[i+1] <- Mw$M1*Mw$K.cold.wet + Mw$M2*Mw$K.warm.arid
    names(Mexp)[i+1] <- Clim.i
  }
  
  Mexp2 <- Mexp[Plot.y]
  Mexp2$K.warm.arid <- Mw$K.warm.arid
  Mexp2$K.cold.wet <- Mw$K.cold.wet
  Mexp2$Pred.cluster <- ifelse(Mw$K.warm.arid >= 0.5, "K-warm/arid", "K-cold/wet")
  Mexp <- list(Mexp, Mexp2)
  
  if(is.null(Save.path) == F){saveRDS(Mexp, Save.path)}
  
  #### Export plot ####
  if(is.null(Save.plot) == F){ggsave(filename = Save.plot, p.main, width = W*0.026458333, height = H*0.026458333, units = "cm")}
  if(return.plot == T){return(p.main)}
  if(return.plot == F){return(Mexp)}
  
}

Residual.plot <- function(MML = NULL, FT = "BRT", Msurf = NULL, MML.surf.train = NULL, MML.clim.train = NULL, Hexa = T, Dot.size = 1, Calcul.RMSE = F, Manual.titre = NULL,
                          Add.RMSE = F, CV.method = "self", Nb.kfold = 10, Res.bottom = T, Res.right = T, Save.residual = NULL, Add.nb = T, Residual.lims = NULL,
                          Nrow.plot = NULL, Annot.position = "bottomright",
                          Mclim = NULL, Annot = NULL, Save.RMSE = NULL, return.plot = T, return.RMSE = F, model = NULL){
  #### Global settings ####
  if(is.null(Msurf) == F & is.null(Mclim) == F){
    if(identical(row.names(Msurf), row.names(Mclim)) == F){
      print("**** Since you want internal CV, we removed the surface sample not used to train the BRT, in the Mclim data.frame. ****")
      Mclim <- Mclim[row.names(Mclim) %in% row.names(Msurf),]}
  }
  
  #### External CV function ####
  External.CV.calc <- function(FT = "BRT", Msurf, Mclim, MML){
    #### Settings ML ####
    if(FT == "BRT"){
      print("****External CV is not available from dismo(), we use gbm() to plot the obs. vs. pred. plot.****")
      training_ml <- function(train_MML.CV, MML, response_var){
        set.seed(123)
        gbm_model <- gbm(
          formula = as.formula(paste(response_var, "~ .")),
          data = train_MML.CV, distribution = MML[[response_var]]$distribution, n.trees = MML[[response_var]]$n.trees,
          interaction.depth = MML[[response_var]]$interaction.depth, shrinkage = MML[[response_var]]$shrinkage, n.minobsinnode = MML[[response_var]]$n.minobsinnode,
          bag.fraction = MML[[response_var]]$bag.fraction, train.fraction = MML[[response_var]]$train.fraction, verbose = F)
        return(gbm_model)
      }
      
      predict_ml <- function(gbm_model, train_MML.CV, MML, response_var){
        preds <- predict(gbm_model, newdata = test_MML.CV, n.trees = MML[[response_var]]$n.trees)
        return(preds)
      }
    }
    
    if(FT == "RF"){
      library(gbm)
      training_ml <- function(train_MML.CV, MML, response_var){
        rf_model <- randomForest(as.formula(paste(response_var, "~ .")), data = train_MML.CV, ntree = MML[[response_var]]$ntree, mtry = MML[[response_var]]$mtry, na.action = na.roughfix, verbose = F)
        return(rf_model)
      }
      
      predict_ml <- function(gbm_model, train_MML.CV, MML, response_var){
        preds <- predict(gbm_model, newdata = test_MML.CV, n.trees = MML[[response_var]]$ntree)
        return(preds)
      }
    }
    
    #### Settings others ####
    library(caret)
    set.seed(123)
    if(is.null(MML.surf.train) == T){MML.surf.train <- Msurf}
    if(identical(MML.surf.train, Msurf) | is.null(Msurf) == T){Combine.intern.extern = F}
    else{Combine.intern.extern = T}
    
    if(is.null(MML.clim.train) == T){
      MML.clim.train <- Mclim
      if(identical(row.names(MML.surf.train), row.names(MML.clim.train)) == F){MML.clim.train <- MML.clim.train[row.names(MML.clim.train) %in% row.names(MML.surf.train),]}
    }
    
    response_var <- names(MML)[1]
    MML.CV <- cbind(MML.surf.train, MML.clim.train[response_var])
    folds <- createFolds(MML.CV[[response_var]], k = Nb.kfold, list = TRUE)
    cv_preds <- rep(NA, nrow(MML.CV))
    
    #### Loop for k-folds ####
    for (i in seq_along(folds)) {
      test_idx <- folds[[i]]
      train_MML.CV <- MML.CV[-test_idx,]
      test_MML.CV <- MML.CV[test_idx,]
      gbm_model <- training_ml(train_MML.CV, MML, response_var) 
      cv_preds[test_idx] <- predict_ml(gbm_model, train_MML.CV, MML, response_var)
    }
    
    #### Combine k-fold (shuffled training/testing set + independant other testing set) ####
    if(Combine.intern.extern == T){
      Msurf_out <- Msurf[!row.names(Msurf) %in% row.names(MML.surf.train),]
      Mclim_out <- Mclim[row.names(Mclim) %in% row.names(Msurf_out),]
      cv_preds_out <- predict(MML[[response_var]], newdata = Msurf_out, n.trees = MML[[response_var]]$n.trees)
      cv_preds <- c(cv_preds, cv_preds_out)
      Mclim <- rbind(MML.clim.train, Mclim_out)
    }
    
    return(list(cv_preds = cv_preds, Mclim = Mclim))
  }
  
  #### Residual calculation ####
  if(is.null(MML) == F){
    #### BRT ####
    if(FT == "BRT"){
      if(CV.method == "external"){
        BRT.res <- External.CV.calc(FT = FT, Msurf, Mclim, MML)
        cv_preds <- BRT.res$cv_preds
        Mclim <- BRT.res$Mclim
        Mconfu1 <- list(BRT = list(cv_preds))
      }
      if(CV.method == "internal"){
        if(is.null(Msurf) == T){Msurf <- MML.surf.train}
        if(is.null(Mclim) == T){Mclim <- MML.clim.train}
        if(nrow(Msurf) == MML$Best.Param$Nb[1]){
          if(is.null(MML[[1]]$fold.fit) == F){Mconfu1 <- list(BRT = list(MML[[1]]$fold.fit))}
          else{
            print("**** Internal CV impossible! The CV result from internal k-fold were not recorded in the BRT list. ****")
            CV.method = "self"
          }
        }
        else{
          print("**** Internal CV impossible! You need to add the same Mpol dataset that the one you used to train the BRT. ****")
          CV.method = "self"
        }
      }
      if(CV.method == "self"){
        if(nrow(Msurf) == MML$Best.Param$Nb[1]){Mconfu1 <- list(BRT = list(MML[[1]]$fit))}
        else{Mconfu1 <- FT.core(Model.BRT = MML, MCore = Msurf, Only.fit = F, Save.RDS = F, Displot = F)
        print("**** Self CV but re-calculating R2 and RMSE from predict(). ****")
        }
      }
      names(Mconfu1$BRT) <- names(MML)[1]
      if(is.null(Msurf) == T){Msurf <- MML.surf.train}
      if(is.null(Mclim) == T){Mclim <- MML.clim.train}
    }
    
    #### RF ####
    if(FT == "RF"){
      if(CV.method == "external"){
        RF.res <- External.CV.calc(FT = FT, Msurf, Mclim, MML) 
        cv_preds <- RF.res$cv_preds
        Mclim <- RF.res$Mclim
        Mconfu1 <- list(RF = list(cv_preds))
      }
      if(CV.method == "internal"){
        if(is.null(Msurf) == T){Msurf <- MML.surf.train}
        if(is.null(Mclim) == T){Mclim <- MML.clim.train}
        if(nrow(Msurf) == MML$Best.Param$Nb[1]){Mconfu1 <- list(RF = list(MML[[1]]$predicted))}
        else{
          print("**** Internal CV impossible! You need to add the same Mpol dataset that the one you used to train the BRT. ****")
          CV.method = "self"
        }
      }
      if(CV.method == "self"){
        Mconfu1 <- FT.core(Model.RF = MML, MCore = Msurf, Only.fit = F, Save.RDS = F, Displot = F)
        print("**** Self CV but re-calculating R2 and RMSE from predict(). ****")
      }
      names(Mconfu1$RF) <- names(MML)[1]
      if(is.null(Msurf) == T){Msurf <- MML.surf.train}
      if(is.null(Mclim) == T){Mclim <- MML.clim.train}
    }
    #### MAT ####
    if(FT == "MAT"){
      Mconfu1 <- FT.core(Model.MAT = MML, MCore = Msurf, Only.fit = F, Save.RDS = F, Displot = F)
    }
    
    #### WAPLS ####
    if(FT == "WAPLS"){
      Mconfu1 <- FT.core(Model.WAPLS = MML, MCore = Msurf, Only.fit = F, Save.RDS = F, Displot = F)
    }
    
    Mconfu1 <- do.call(rbind, lapply(Mconfu1, as.data.frame))
    Mconfu1$Sites <- row.names(Mconfu1)
    Mconfu1 <- setNames(data.frame(reshape2::melt(Mconfu1, id = "Sites")), c("Sites", "Pred.clim.lab", "Pred.clim.val"))
    
    Mconfu1.clim <- Mclim[match(levels(Mconfu1$Pred.clim.lab), names(Mclim))]
  }
  else{
    Mconfu1 <- Msurf
    Mconfu1$Sites <- row.names(Msurf)
    Mconfu1 <- setNames(data.frame(reshape2::melt(Mconfu1, id = "Sites")), c("Sites", "Pred.clim.lab", "Pred.clim.val"))
    Mconfu1$Calib <- Mconfu1$Pred.clim.lab
    Mconfu1$Pred.clim.lab <- gsub("_.*", "", Mconfu1$Pred.clim.lab)
    
    Mconfu1.clim <- Mclim[match(unique(Mconfu1$Pred.clim.lab), names(Mclim))]
  }
  Mconfu1.clim$Sites <- unique(Mconfu1$Sites)
  Mconfu1.clim <- setNames(data.frame(reshape2::melt(Mconfu1.clim, id = "Sites")), c("Sites", "Obs.clim.lab", "Obs.clim.val"))
  
  Mconfu1 <- cbind(Mconfu1, Mconfu1.clim[3])
  Mconfu1$Residuals <- (Mconfu1$Obs.clim.val - Mconfu1$Pred.clim.val)
  
  if(is.null(MML) == T){
    Mconfu1 <- subset(Mconfu1, select = -c(Pred.clim.lab))
    names(Mconfu1)[names(Mconfu1) == "Calib"] <- "Pred.clim.lab"
    Mconfu1 <- Mconfu1[c(1,3,2,4,5)]
  }
  
  #### Residual scaling ####
  Mconfu1_scaled <- Mconfu1 %>%
    dplyr::group_by(Pred.clim.lab) %>%
    mutate(Residuals.sc = (Residuals - mean(Residuals)) / sd(Residuals))
  
  #### RMSE calculation ####
  if(Calcul.RMSE == T | Add.RMSE == T){
    MRes <- Mconfu1[c("Sites", "Pred.clim.lab", "Residuals")]
    MRes <- reshape2::dcast(MRes, Sites ~ Pred.clim.lab, value.var = "Residuals")
    MRes <- sqrt(colSums(subset(MRes, select = - c(Sites))^2)/nrow(MRes))
    
    if(is.null(Save.RMSE) == F){
      print("Local RMSE for each calibration:")
      print(MRes)
      saveRDS(MRes, Save.RMSE)
    }
  }
  else{MRes <- NULL}
  
  #### Lim settlement ####
  Mlims <- Mconfu1[names(Mconfu1) %in% c("Pred.clim.lab", "Pred.clim.val", "Obs.clim.val")]
  Mlims <- rbind(setNames(Mlims[c(1,2)], c("Pred.clim.lab", "B")), setNames(Mlims[c(1,3)], c("Pred.clim.lab", "B")))
  Mlims <- Mlims %>%
    dplyr::group_by(Pred.clim.lab) %>%
    dplyr::summarize(
      xlim = min(B),
      ylim = max(B)
    )
  Mlims <- reshape2::melt(Mlims, id = "Pred.clim.lab")
  Mlims$Sites <- "Dummy.dot"
  Mlims <- subset(Mlims, select = -c(variable))
  Mlims$Obs.clim.val <- Mlims$value
  Mlims$Obs.clim.lab <- Mlims$Pred.clim.lab
  Mlims$Residuals <- 0
  names(Mlims)[names(Mlims) == "value"] <- "Pred.clim.val"
  
  Mlims <- Mlims[match(names(Mconfu1), names(Mlims))]
  Mlims$Residuals.sc <- 0
  
  #### Global graphical param ####
  New.lab1 <- New.lab[match(levels(Mconfu1$Pred.clim.lab), names(New.lab))]
  
  if(all(is.na(unique(New.lab1)))){
    New.lab1 <- levels(Mconfu1$Pred.clim.lab)
    New.lab1 <- paste(sub("_", "~(", New.lab1), ")", sep = "")
    New.lab1 <- sub("_", "~", New.lab1)
    New.lab1 <- gsub("5Me", "[5*Me]", New.lab1)
    New.lab1 <- gsub("~\\[", "[", New.lab1)
  }
  
  if(is.null(model) == F){
    Ytitle <- paste(model, "predictions")
  }
  else{Ytitle <- "Predicted \nclimate parameters"}
  
  if(Hexa == T){My_dots <- geom_hex(bins = 30, aes(fill = ..count..))}
  else{My_dots <- geom_point(color = "grey50", size = Dot.size)}
  
  New.lab <- c("d13Corg" = "paste(delta^13,C[TOC])", "d13Ctot" = "paste(delta^13,C[tot])", "d15Ntot" = "paste(delta^15,N[tot])", "A: Core sections" = "MasterCore",
               "Corg.N_atom" = "C/N", "X.C_inorg_DeltaM" = "TIC", "wt..Ntot" = "wt~'%N'[tot]",
               "Fe/Al" = "Fe/Al", "Rb/Sr" = "Rb/Sr", "Ti/Al" = "Ti/Al", "K/Ca" = "K/Ca", "C14" = "paste(Dating^14,C)",
               "L" = "L\\*", "RABD660_670" = "RABD[660-670]", "Q7.4" = "Q[700/400]", "Chlo_a" = "Chlo[a]~(mg.g^-1)",
               "MAAT" = "MAAT (degree*C)", "MAF" = "MAF (degree*C)", "MAP" = "MAP~(mm.yr^-1)",  "MPCOQ" = "MPCOQ~(mm.yr^-1)", "Altitude" = "Altitude", "AI" = "AI", "Pspr" = "P[spring]~(mm.yr^-1)",
               "IIIa.IIa" = "Sigma(IIIa/IIa)", "Ib.Ia" = "Ib/Ia", "IR6_7Me" = "IR[6+7~Me]","IRp6_7Me" = "IRp[6+7~Me]", "GDGT0.Crenar" = "GDGT[0]/Crenar", "pCren" = "\'%\'[Cren]",
               "Art_Ama/Poa" = "over((Artemisia+Amaranth.), Poaceae)", "Poa/Ama" = "over(Poaceae, Amaranthaceae)", "Ama/Art" = "over(Amaranthaceae,Artemisia)", "Art/Ama+Art" = "over(Artemisia,(Amaranthaceae+Artemisia))", "Art+Ama/Poa" = "over((Artemisia+Amaranthaceae),Poaceae)", "Poa/Art" = "over(Poaceae,Artemisia)", "AP/NAP" = "over(Arboreal~Pollen, Non-Arboreal~Pollen)",
               "PC1" = "PC1[XRF]", "PC2" = "PC2[XRF]", "MS" = "MS~(10^-5~SI)", "Sus_mag" = "")
  
  if(is.null(Residual.lims) == F){
    Res.lim.r <- xlim(Residual.lims)
    Res.lim.b <- ylim(Residual.lims)
    
  }
  else{Res.lim.b <- NULL; Res.lim.r <- NULL}
  
  if(Res.bottom == T){Myaxis.x <- element_blank()}
  else{Myaxis.x <- element_line()}
  
  #### Plotting loops ####
  plots <- list()
  for(i in 1:length(unique(Mconfu1$Pred.clim.lab))){
    #### Param i ####
    Param <- unique(Mconfu1$Pred.clim.lab)[i]
    New.lab1.i <- New.lab1[i]
    if(is.null(Annot) == F){New.lab1.i <- paste("(", Annot[i], ")~", New.lab1.i, sep = "")}
    Mconfu.i <- Mconfu1[Mconfu1$Pred.clim.lab == Param,]
    Mconfu1_scaled.i <- Mconfu1_scaled[Mconfu1_scaled$Pred.clim.lab == Param,]
    
    if(i == 1){Ylab <- element_text()}
    if(i > 1){Ylab <- element_blank()}
    
    #### Add Annotation ####
    if(Add.RMSE == T){
      R2 <- summary(lm(Pred.clim.val ~ Obs.clim.val, data = Mconfu.i))$r.squared
      
      label_text <- paste0("R^2==", signif(R2, 2), "*','~RMSE==", signif(MRes[[i]], 2))
      
      if(Annot.position == "bottomright"){X1 = Inf; Y1 = -Inf; H1 = 1.1; V1 = -1}
      if(Annot.position == "bottomleft"){X1 = -Inf; Y1 = -Inf; H1 = 1.1; V1 = -1}
      if(Annot.position == "topright"){X1 = Inf; Y1 = Inf; H1 = 1.1; V1 = -1}
      if(Annot.position == "topleft"){X1 = -Inf; Y1 = Inf; H1 = -0.1; V1 = 1.5}
      
      if(Add.nb == T){
        n_obs <- length(Mconfu.i$Obs.clim.val)
        label_text <- paste0("R^2==", signif(R2, 2), "*','~RMSE==", signif(MRes[[i]], 2), "*','~n==", n_obs)
      }
      Annotate.R2 <- annotate("text", x = X1, y = Y1, label = label_text, hjust = H1, vjust = V1, parse = TRUE, size = 3)
    }
    else{Annotate.R2 <- stat_poly_eq(size = 3, vstep = 0.07, label.y = "bottom", label.x = "right")}
    
    #### Manual title ####
    if(is.null(Manual.titre) == F){
      if(is.null(Annot) == F){Manual.titre <- paste("(", Annot[i], ")~", Manual.titre, sep = "")}
      
      New.lab1.i[1] <- Manual.titre
    }
    
    #### Plot (main) ####
    p.main <- ggplot(Mconfu.i, aes(x = Obs.clim.val, y = Pred.clim.val)) +
      My_dots +
      geom_point(data = Mlims, color = NA, fill = NA)+
      geom_smooth(method = "lm", se = F, linewidth = 0.7, formula = 'y ~ x') +
      scale_fill_gradientn(name = "Nb. obs.",
                           colors = c("grey90", "royalblue", "darkorange", "darkred"),  # Discrete color steps
                           values = scales::rescale(c(0, 1, 10, 100)),  # Scale the breaks
                           breaks = c(0, 1, 10, 100)) +
      geom_abline(slope = 1, intercept = 0, color = "grey20", linetype = "dashed") + 
      labs(x = NULL, 
           y = Ytitle) + 
      ggtitle(parse(text = New.lab1.i[1]))+
      Annotate.R2 +
      theme(panel.background = element_rect(fill = NA, colour = 'black', linewidth = .5), panel.border = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 13),
            legend.position = "none", panel.grid = element_line(colour = "grey80", linetype = 2, linewidth = .1),
            axis.text.x = element_blank(),
            axis.title.y = Ylab, axis.ticks.x = Myaxis.x,
            axis.line.x = element_blank())
    
    #### Plot bottom (residuals) ####
    if(Res.bottom == T){
      presidu.b <- ggplot(Mconfu1_scaled.i, aes(x = Obs.clim.val, y = Residuals)) +
        geom_point(aes(color = Residuals), size = Dot.size) +
        geom_point(data = Mlims, color = NA, fill = NA)+
        geom_hline(yintercept = 0, lty = "dashed")+
        labs(x = "Observed climate parameters",
             y = "Residuals") + Res.lim.b +
        scale_color_gradient2(name = "Residual (z-scores)",
                              low = "#963327",
                              mid = "white",
                              high = "#963327",
                              midpoint = 0,
                              limit = c(-max(abs(Mconfu1_scaled.i$Residuals)), max(abs(Mconfu1_scaled.i$Residuals))))+
        theme(
          plot.background = element_blank(), plot.margin = unit(c(0,0,0,0), 'pt'),
          axis.title.y = Ylab,
          legend.position = "none", panel.grid = element_line(colour = "grey80", linetype = 2, linewidth = .1),
          panel.background = element_rect(fill = NA, colour = 'black'), panel.border = element_blank(),
          strip.text = element_blank(), strip.background = element_blank())
      
      p <- (p.main / presidu.b)+ plot_layout(nrow = 2, heights = c(.8,.2))
    }
    
    #### Plot right (residuals) ####
    if(Res.right == T){
      presidu.r <- ggplot(Mconfu1_scaled.i, aes(y = Pred.clim.val, x = Residuals)) +
        geom_point(aes(color = Residuals), size = Dot.size) +
        geom_point(data = Mlims, color = NA, fill = NA)+
        geom_vline(xintercept = 0, lty = "dashed")+
        labs(y = "Predicted climate parameters",
             x = "Residuals") + Res.lim.r +
        scale_color_gradient2(name = "Residual (z-scores)",
                              low = "#963327",
                              mid = "white",
                              high = "#963327",
                              midpoint = 0,
                              limit = c(-max(abs(Mconfu1_scaled.i$Residuals)), max(abs(Mconfu1_scaled.i$Residuals))))+
        theme(
          plot.background = element_blank(),  plot.margin = unit(c(0,0,0,0), 'pt'),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "none", panel.grid = element_line(colour = "grey80", linetype = 2, linewidth = .1),
          panel.background = element_rect(fill = NA, colour = 'black'), panel.border = element_blank(),
          strip.text = element_blank(), strip.background = element_blank())
      
      if(Res.bottom == F){p <- (p.main | presidu.r)+ plot_layout(ncol = 2, widths = c(.75,.25))}
      else{
        p.dummy <- ggplot(Mconfu1_scaled.i, aes(y = Pred.clim.val, x = Residuals)) + Res.lim.b + Res.lim.r +
          theme(plot.background = element_blank(),  plot.margin = unit(c(0,0,0,0), 'pt'), axis.title = element_blank(), panel.border = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.grid = element_blank(),
                panel.background = element_blank(), strip.text = element_blank(), strip.background = element_blank())
        
        p.main <- p.main + theme(plot.margin = ggplot2::margin(0,0,0,0))
        presidu.r <- presidu.r + theme(plot.margin = ggplot2::margin(0,0,0,0))
        presidu.b <- presidu.b + theme(plot.margin = ggplot2::margin(0,0,0,0))
        p.dummy <- p.dummy + theme(plot.margin = ggplot2::margin(0,0,0,0))
        
        p <- ((p.main | presidu.r) + plot_layout(widths = c(0.8, 0.2))) /
          ((presidu.b | p.dummy) + plot_layout(widths = c(0.8, 0.2))) +
          plot_layout(heights = c(0.8, 0.2)) # only heights at top level
      }
    }
    
    if(Res.bottom == F & Res.right == F){p <- p.main}
    plots[[i]] <- p
  }
  
  if(is.null(Nrow.plot) == F){plots <- wrap_plots(plots, nrow = Nrow.plot)}
  else{plots <- wrap_plots(plots)}
  
  #### Export ####
  if(is.null(Save.residual) == F){saveRDS(Mconfu.i, Save.residual)}
  if(return.plot == T){return(plots)}
  if(return.RMSE == T){return(MRes)}
  if(return.RMSE == F && return.plot == F){return(Mconfu.i)}
}

RMSE.comb <- function(ML.comb, BRT.K1, BRT.K2) {
  Average.wg <- colMeans(ML.comb[[2]][2:3])
  RMSE.KAR <- BRT.K1$Best.Param$RMSE.choice
  RMSE.KWT <- BRT.K2$Best.Param$RMSE.choice
  A <- data.frame(sqrt(Average.wg[1]^2*RMSE.KAR^2 + Average.wg[2]^2*RMSE.KWT^2))
  row.names(A) <- row.names(BRT.K1$Best.Param)
  names(A) <- "RMSE.choice"
  A <- list(Best.Param = A)
  return(A)
}

GDGT.ensemble.RMSE <- function(Msurf, Param.clim, Ensemble, Mclim, Mpaleo, Mpaleo.max, Mpaleo.min, Variance.method = "quadratic",
                               Ensemble.WM = T, Ensemble.BAY = T, Save.path = NULL){
  
  for (i in 1:length(Param.clim)){
    #### RMSE calculation ####
    Param.clim.i <- Param.clim[i]
    Ensemble.i <- Ensemble[[i]]
    M.RMSE.surf <- Msurf[grep(paste("^", Param.clim.i, sep = ""), names(Msurf))]
    PB <- setdiff(Ensemble.i, names(M.RMSE.surf))
    if(length(PB) > 0){
      print(paste("The following model is missing in the surface data:", PB))
      next
    }
    RMSE <- Residual.plot(MML = NULL, Msurf = M.RMSE.surf[Ensemble.i], Mclim = Mclim, return.plot = F, Calcul.RMSE = T, return.RMSE = T)
    
    #### Weighted variance function ####
    row_weighted_variance <- function(mat, weights, method = c("standard", "unbiased", "quadratic")) {
      method <- match.arg(method)
      
      if (ncol(mat) != length(weights)) {
        stop("Number of columns in matrix must match length of weights.")
      }
      
      apply(mat, 1, function(x) {
        w_mean <- sum(weights * x) / sum(weights)
        
        if (method == "standard") {
          return(sum(weights * (x - w_mean)^2) / sum(weights))
          
        } else if (method == "unbiased") {
          w_sum <- sum(weights)
          w_squared_sum <- sum(weights^2)
          correction = w_sum - (w_squared_sum / w_sum)
          return(sum(weights * (x - w_mean)^2) / correction)
          
        } else if (method == "quadratic") {
          return(sum(weights^2 * (x - w_mean)^2))
        }
      })
    }
    
    #### Weighted Mean RMSE ####
    if(Ensemble.WM == T){
      wg <- 1/RMSE
      wg <- wg/sum(wg)
      Mpaleo[paste(Param.clim.i, "Ensemble", sep = "_")] <- as.numeric(as.matrix(Mpaleo[Ensemble.i]) %*% wg)
      prediction_matrix <- as.matrix(Mpaleo[Ensemble.i])
      weighted_variance <-row_weighted_variance(prediction_matrix, wg, method = Variance.method)
      
      Mpaleo.max[paste(Param.clim.i, "Ensemble", sep = "_")] <- Mpaleo[paste(Param.clim.i, "Ensemble", sep = "_")] + 1.96*sqrt(weighted_variance)
      Mpaleo.min[paste(Param.clim.i, "Ensemble", sep = "_")] <- Mpaleo[paste(Param.clim.i, "Ensemble", sep = "_")] - 1.96*sqrt(weighted_variance)
    }
    
    #### Bayesian RMSE ####
    if(Ensemble.BAY == T){
      # RMSE to approximate log-likelihoods (assuming Gaussian error)
      # log_likelihood <- -nrow(Mpaleo)/2 * log(RMSE.MAAT^2)  # Simplified proxy; more precise would be 
      log_likelihood <- -1/2 * log(RMSE^2)  # Simplified proxy; more precise would be 
      
      # Convert to posterior model probabilities (softmax)
      posterior_weights <- exp(log_likelihood - max(log_likelihood))  # for numerical stability
      posterior_weights <- posterior_weights / sum(posterior_weights)
      Mpaleo[paste(Param.clim.i, "Ensemble_BAY", sep = "_")] <- as.numeric(as.matrix(Mpaleo[Ensemble.i]) %*% posterior_weights)
      
      prediction_matrix <- as.matrix(Mpaleo[Ensemble.i])
      posterior_variance <- apply(prediction_matrix, 1, function(x) {sum(posterior_weights * (x - sum(posterior_weights * x))^2)})
      
      se <- sqrt(posterior_variance)
      Mpaleo.max[paste(Param.clim.i, "Ensemble_BAY", sep = "_")] <- Mpaleo[paste(Param.clim.i, "Ensemble_BAY", sep = "_")] + 1.96*sqrt(posterior_variance)
      Mpaleo.min[paste(Param.clim.i, "Ensemble_BAY", sep = "_")] <- Mpaleo[paste(Param.clim.i, "Ensemble_BAY", sep = "_")] - 1.96*sqrt(posterior_variance)}
    
  }
  #### Export ####
  if(is.null(Save.path) == F){
    saveRDS(Mpaleo, Save.path)
    saveRDS(Mpaleo.max, gsub("\\.Rds", "_MaxI.Rds", Save.path))
    saveRDS(Mpaleo.min, gsub("\\.Rds", "_MinI.Rds", Save.path))
  }
  
  Mexp <- list(Mpaleo)
  if(is.null(Mpaleo.min) == F){Mexp[[length(Mexp) + 1]] <- Mpaleo.min}
  if(is.null(Mpaleo.max) == F){Mexp[[length(Mexp) + 1]] <- Mpaleo.max}
  return(Mexp)
}
