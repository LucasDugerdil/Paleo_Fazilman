#### Label settings ####
New.lab <- c("d13Corg" = "paste(delta^13,C[TOC])", "d13Ctot" = "paste(delta^13,C[tot])", "d15Ntot" = "paste(delta^15,N[tot])", "A: Core sections" = "MasterCore",
             "Corg.N_atom" = "C/N", "X.C_inorg_DeltaM" = "TIC", "wt..Ntot" = "wt~'%N'[tot]",
             "Fe/Al" = "Fe/Al", "Rb/Sr" = "Rb/Sr", "Ti/Al" = "Ti/Al", "K/Ca" = "K/Ca", "C14" = "paste(Dating^14,C)",
             "L" = "L\\*", "RABD660_670" = "RABD[660-670]", "Q7.4" = "Q[7/4]", "Chlo_a" = "Chlo[a]~(mg.g^-1)",
             "MAAT" = "MAAT", "MAP" = "MAP~(mm.yr^-1)", "Altitude" = "Altitude", "AI" = "AI", "Pspr" = "P[spring]~(mm.yr^-1)",
             "IIIa.IIa" = "Sigma(IIIa/IIa)", "Ib.Ia" = "Ib/Ia", "IR6_7Me" = "IR[6+7~Me]","IRp6_7Me" = "IRp[6+7~Me]","IR6Me" = "IR[6~Me]","IR7Me" = "IR[7~Me]", "GDGT0.Crenar" = "GDGT[0]/Crenar", "pCren" = "\'%\'[Cren]",
             "pCren.p" = "\'%\'[Cren*minute]",
             # "Art_Ama/Poa" = "(Ar.+Amar.)/Poa.", "Poa/Ama" = "Poa./Amar.", "Ama/Art" = "Amar./Ar.", "Art/Ama+Art" = "Ar./(Amar.+Ar.)", "Art+Ama/Poa" = "(Ar.+Amar.)/Poa.", "Poa/Art" = "Poa./Ar.",
             "Art_Ama/Poa" = "over((Artemisia+Amaranth.), Poaceae)", "Poa/Ama" = "over(Poaceae, Amaranthaceae)", "Ama/Art" = "over(Amaranthaceae,Artemisia)", "Art/Ama+Art" = "over(Artemisia,(Amaranthaceae+Artemisia))", "Art+Ama/Poa" = "over((Artemisia+Amaranthaceae),Poaceae)", "Poa/Art" = "over(Poaceae,Artemisia)", "AP/NAP" = "over(Arboreal~Pollen, Non-Arboreal~Pollen)",
             "PC1" = "PC1[XRF]", "PC2" = "PC2[XRF]", "MS" = "MS~(10^-5~SI)", "Sus_mag" = "")

My.unit <- c("K/Ca" = NA, "Rb/Sr" = NA, "Rb/K" = NA, "Fe/Mn" = NA, "Ca/Mg" = NA, "Ca/Ti" = NA, "Al/Si" = NA, "Zr/Ti" = NA, "A: Core sections" = NA, "Inc/Coh" = NA, "Mn/Fe" = NA, "Coh/Inc" = NA, "Fe/Ti" = NA,
             "Ti/Al" = NA, "Fe/Al" = NA, "Br/Ti" = NA, "Sr/Ca" = NA, "S/Fe" = NA, "Si/Ti" = NA, "S/Ti" = NA, "K/Ti" = NA, "Mg/Ca" = NA, "Sr/Ti" = NA, "Coh" = NA, "Inc" = NA, "C14" = NA, "MI" = NA,
             "Al" = NA, "Si" = NA, "K" = NA, "Ti" = NA, "S" = NA, "Mn" = NA, "Fe" = NA, "Ca" = NA, "Zr" = NA, "Zn" = NA, "Rb" = NA, "Sr" = NA, "Br" = NA, "Cu" = NA, "Pb" = NA, "Ar" = NA, "P" = NA,
             "PC1" = NA, "PC2" = NA, "BIT" = "%", "IIIa.IIa" = NA, "Q7.4" = NA, "RABD660_670" = NA, "CONISS" = NA, "Sus_mag" = "10e-5 SI", "Chlo_a" = NA, "Density" = NA, "MS" = NA, 
             "Ib.Ia" = NA, "pCren" = NA, "pCren.p" = NA, "IR6_7Me" = "%", "IRp6_7Me" = "%", "GDGT0.Crenar" = NA, "IR" = NA, "CI" = NA, " " = NA,
             "MAAT" = "°C", "MAP" = NA, "Altitude" = "m a.s.l.", "AI" = NA, "Pspr" = NA, "Units" = NA,
             "Poa/Art" = NA, "Poa/Ama" = NA, "Ama/Art" = NA, "Art/Ama+Art" = NA, "AP/NAP" = NA, "Art_Ama/Poa" = NA,
             "Corg.N_atom" = NA, "d13Corg" = "%", "d15Ntot" = "%", "d13Ctot" = "%")

Fazilman.clus.color <- c("#7042b9ff", "#6872b3ff", "#35978F", "#737591ff", "grey75", "#c1982eff", "#BF812D", "#4b3f0dff")

TabRa <- list(data.frame(Deno = c("Artemisia"), Nomi = c("Chenopodiaceae")),
              data.frame(Deno = c("Artemisia", "Chenopodiaceae"), Nomi = c("Poaceae"))
              )

#### Prentice colors ####
Couleur.Prentice <- c(
  "WAMX"="#185699FF",
  "TEDE"="#72a72bff",
  "XERO"="#cf0156ff",
  "COMX"="#c1b646ff",
  "HODE"="#fcf8b6ff",
  "WAST"="#e2a064ff",
  "CLMX"="#ee9b2fff",
  "PION"="#0e3056ff",
  "TAIG"="#32156eff",
  "TUND"="#caa6c2ff",
  "COCO"="#9d2e58ff",
  "COST"="#f3c768ff",
  "ANTH"="#6e1d2cff",
  "CODE"="#d6d81eff",
  "CLDE"="#b1d5f0ff",
  "No_data" = "grey80")

Couleur.Prentice <- Couleur.Prentice[order(names(Couleur.Prentice))]

#### Functions ####
Plot.geoch.core <- function(XRF, MS, Isotope, GDGT, Select.interv.x, XRF.SD, GDGT.min, GDGT.max, 
                            Hiatus.age, Axis.x.tible, Nb.ticks.minors, Change.titles = NULL,
                            Plot.x, CONISS, Nzone, Ratio, Limites.x, Remove.GDGT, Keep.GDGT, 
                            Print.cluster.zone, Reverse.arrow, Smooth.sd = T, Smooth.proxy = NULL,
                            Remove.MS, Groupes, Remove.Isotope, Keep.Isotope, Manual.line, return.plot = F, 
                            Display.legend, Arrows, Keep.MS, Nsous.zone, Manual.small.line, Manual.color.curves = NULL,
                            Col.area, Cparam, Remove.XRF, Keep.XRF, Save.CONISS, Smooth.param = 0.6,
                            Zone.clim2, Name.zone2, Temp.zone2, Panel.param = NULL, Arrow.y = NULL,
                            Save.plot, W, H, Zone.clim, Name.zone, Temp.zone, C14.path = NULL){
  #### Initialisation ####
  library(tidypaleo)
  library(ggrepel)
  if(missing(Nb.ticks.minors)){Nb.ticks.minors = NULL}
  if(missing(Axis.x.tible)){Axis.x.tible = NULL}
  if(missing(Hiatus.age)){Hiatus.age = NULL}
  if(missing(XRF)){XRF = NULL}
  if(missing(XRF.SD)){XRF.SD = NULL}
  if(missing(Isotope)){Isotope = NULL}
  if(missing(MS)){MS = NULL}
  if(missing(GDGT)){GDGT = NULL}
  if(missing(GDGT.min)){GDGT.min = NULL}
  if(missing(GDGT.max)){GDGT.max = NULL}
  if(missing(Plot.x)){stop("Select the age / depth feature for the plot.", call. = F)}
  if(missing(CONISS)){CONISS = F}
  if(missing(Col.area)){Col.area = NULL}
  if(missing(Arrows)){Arrows = F}
  if(missing(Display.legend)){Display.legend = F}
  if(missing(Print.cluster.zone)){Print.cluster.zone = F}
  if(missing(Manual.small.line)){Manual.small.line = NULL}
  if(missing(Ratio)){Ratio = NULL}
  if(missing(Groupes)){Groupes = NULL}
  if(missing(Nzone)){Nzone = NULL}
  if(missing(Nsous.zone)){Nsous.zone = NULL}
  if(missing(Manual.line)){Manual.line = NULL}
  if(missing(Remove.XRF)){Remove.XRF = NULL}
  if(missing(Remove.Isotope)){Remove.Isotope = NULL}
  if(missing(Remove.MS)){Remove.MS = NULL}
  if(missing(Remove.GDGT)){Remove.GDGT = NULL}
  if(missing(Keep.XRF)){Keep.XRF = NULL}
  if(missing(Keep.MS)){Keep.MS = NULL}
  if(missing(Keep.Isotope)){Keep.Isotope = NULL}
  if(missing(Keep.GDGT)){Keep.GDGT = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Save.CONISS)){Save.CONISS = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Select.interv.x)){Select.interv.x = 1000}
  if(missing(Reverse.arrow)){Reverse.arrow = NULL}
  if(missing(Zone.clim)){Zone.clim = NULL}
  if(missing(Zone.clim2)){Zone.clim2 = NULL}
  if(missing(Name.zone2)){Name.zone2 = NULL}
  if(missing(Name.zone)){Name.zone = NULL}
  if(missing(Temp.zone) & is.null(Zone.clim) == F){Temp.zone = rep("C", length(Zone.clim))}
  if(missing(Temp.zone) & is.null(Zone.clim) == T){Temp.zone = F}
  if(missing(Temp.zone2)){Temp.zone2 = NULL}
  
  Possible.x <- c("Age", "mean", "Litho", "Units", "MCD", "X", "SECT.NUM", "SECT.DEPTH", "MCD.1",
                  "Ordin1", "Ordin2", "Ordin3", "Ordin4", "Ordin5")
  Keep.zone.clim.full <- NULL
  
  #### Color ramp setting for zone temp ####
  if(any(unique(grepl("#", Temp.zone))) == T){
    print("Manual color scale for climate zone activated.")
    Rect.color.scale <- Temp.zone
    names(Rect.color.scale) <- Rect.color.scale
    if(any(unique(grepl("#", Temp.zone2))) == T){
      Rect.color.scale <- unique(c(Temp.zone, Temp.zone2))
      names(Rect.color.scale) <- Rect.color.scale
      
    }
    
  }
  else{
    Rect.color.scale <- c('grey')
    if(length(unique(Temp.zone)) == 0){Rect.color.scale <- c("grey")}
    if(length(unique(Temp.zone)) == 3){Rect.color.scale <- c("#75AADB", "black", "#E76D51")}
    if(length(unique(Temp.zone)) == 1 & "G" %in% unique (Temp.zone)){Rect.color.scale <- c("grey")}
    if(length(unique(Temp.zone)) == 1 & "C" %in% unique(Temp.zone)){Rect.color.scale <- c("#75AADB")}
    if(length(unique(Temp.zone)) == 1 & "W" %in% unique(Temp.zone)){Rect.color.scale <- c("#E76D51")}
    if(length(unique(Temp.zone)) == 2 & "C" %in% unique(Temp.zone) & "W" %in% unique(Temp.zone)){Rect.color.scale <- c("#E76D51", "#75AADB")}
    if(length(unique(Temp.zone)) == 2 & "G" %in% unique(Temp.zone) & "W" %in% unique(Temp.zone)){Rect.color.scale <- c("black", "#E76D51")}
    if(length(unique(Temp.zone)) == 2 & "G" %in% unique(Temp.zone) & "C" %in% unique(Temp.zone)){Rect.color.scale <- c("#75AADB", "black")}
  }
  
  #### Zone temp settings ####
  if(is.null(Zone.clim) == F){
    Keep.zone.clim.full <- Zone.clim
    if(missing(Limites.x) == F){
      Zone.clim <- c(Limites.x[1], Zone.clim[which(Zone.clim >= Limites.x[1] & Zone.clim <= Limites.x[2])], Limites.x[2])
    }
  }
  
  if(is.null(Col.area) == F){
    if((is.null(Zone.clim) == F | is.null(Zone.clim2) == F)){
      if(Col.area == 1){
        My_name_zone = Name.zone
        My_temp_zone = Temp.zone
        My_zone_clim = Zone.clim
      }
      if(Col.area == 2){
        My_name_zone = Name.zone2
        My_temp_zone = Temp.zone2
        My_zone_clim = Zone.clim2
      }
      
      
      yo = data.frame(xmin = My_zone_clim[seq(1,length(My_zone_clim), by=2)], 
                      xmax = My_zone_clim[seq(2,length(My_zone_clim), by=2)], 
                      Temp.col = My_temp_zone)}
    else{yo = data.frame(xmin = 0, xmax = 0, Temp.col = "black"); My_zone_clim <- NULL}
    
    if((is.null(Name.zone) == F | is.null(Name.zone2) == F) & is.null(My_zone_clim) == F){
      yo2 = data.frame(xmin = My_zone_clim[seq(1,length(My_zone_clim), by=2)], 
                       xmax = My_zone_clim[seq(2,length(My_zone_clim), by=2)],
                       Temp.col = My_temp_zone,
                       Title.zone = My_name_zone)}
    else{yo2 = data.frame(xmin = 0, xmax = 0, Temp.col = "", Title.zone = "", Categorie = "")}
    
    Zone.clim.to.plot <- geom_rect(data = yo, inherit.aes = F, na.rm = T,
                                   mapping = aes(xmin = -Inf, xmax= +Inf, ymin = xmin, ymax = xmax, fill = Temp.col),
                                   alpha=0.1, color = "grey", linewidth = 0.2, linetype = 0)
    
  }
  else{Zone.clim.to.plot <- NULL}
  
  if(is.null(Hiatus.age) == F){
    Hiatus.area <- data.frame(xmin = Hiatus.age[1], xmax = Hiatus.age[2])
    Hiatus.area <- geom_rect(data = Hiatus.area, inherit.aes = F,
                             mapping = aes(ymin=xmin, ymax=xmax, xmin=-Inf, xmax=+Inf), fill = "grey97",
                             alpha = 1, color = "grey60", linewidth = 0, linetype = 0, na.rm = T)}
  else{Hiatus.area <- NULL}
  
  #### XRF ####
  if(is.null(XRF) == F){
    XRF.keep <- XRF[names(XRF) %in% Possible.x]
    XRF <- XRF[!names(XRF) %in% setdiff(Possible.x, Plot.x)]
    #### XRF interval ####
    if(is.null(XRF.SD) == F){
      Plot.SD <- Plot.x
      XRF.SD[setdiff(names(XRF), names(XRF.SD))] <- NA
      XRF.SD <- XRF.SD[names(XRF)]
      XRF.max <- cbind(XRF[Plot.SD], XRF[setdiff(names(XRF.SD),Plot.SD)] + 1.96*XRF.SD[setdiff(names(XRF.SD),Plot.SD)])
      XRF.min <- cbind(XRF[Plot.SD], XRF[setdiff(names(XRF.SD),Plot.SD)] - 1.96*XRF.SD[setdiff(names(XRF.SD),Plot.SD)])
      Mratio.max <- XRF[0]
      Mratio.min <- XRF[0]
    }
    
    #### Calcul des Ratios ####
    if(is.null(Ratio) == F){
      Mratio <- XRF[0]
      deno <- gsub("/.*", "", Ratio)
      nomi <- gsub(".*/", "", Ratio)
      for(i in 1:length(Ratio)){
        Mratio[i] <- XRF[deno[i]] / XRF[nomi[i]] 
        names(Mratio)[i] <- Ratio[i]
        
        if(is.null(XRF.SD) == F){
          Mratio.min[i] <- XRF.min[deno[i]] / XRF.min[nomi[i]] 
          Mratio.max[i] <- XRF.max[deno[i]] / XRF.max[nomi[i]] 
          names(Mratio.min)[i] <- Ratio[i]
          names(Mratio.max)[i] <- Ratio[i]
        }
      }
    }
    #### Selection des proxy ####
    if(is.null(Remove.XRF) == F){
      XRF <- XRF[setdiff(names(XRF), Remove.XRF)]
      if(is.null(XRF.SD) == F){
        XRF.min <- XRF.min[setdiff(names(XRF.min), Remove.XRF)]
        XRF.max <- XRF.max[setdiff(names(XRF.max), Remove.XRF)]
      }
    }
    if(is.null(Keep.XRF) == F){
      XRF <- XRF[c(Keep.XRF, Plot.x)]
      if(is.null(XRF.SD) == F){
        XRF.min <- XRF.min[c(Plot.x, intersect(Keep.XRF, names(XRF.min)))]
        XRF.max <- XRF.max[c(Plot.x, intersect(Keep.XRF, names(XRF.max)))]
      }
    }
    
    #### Mise en forme de la matrice ####
    if(is.null(Ratio) == F){
      XRF <- cbind(XRF, Mratio)
      if(is.null(XRF.SD) == F){
        XRF.min <- cbind(XRF.min, Mratio.min)
        XRF.max <- cbind(XRF.max, Mratio.max)
      }
    }
    n.XRF <- ncol(XRF)-1
    
    XRF <- reshape2::melt(XRF, id = Plot.x)
    names(XRF)[1] <- "Plot.x"
    XRF <- na.omit(XRF)                 # iniore les valeurs maquantes
    
    if(is.null(XRF.SD) == F){
      XRF.min <- reshape2::melt(XRF.min, id = Plot.x)
      names(XRF.min)[1] <- "Plot.x"
      names(XRF.min)[3] <- "Min"
      XRF.min <- na.omit(XRF.min)                 # iniore les valeurs maquantes
      
      XRF.max <- reshape2::melt(XRF.max, id = Plot.x)
      names(XRF.max)[1] <- "Plot.x"
      names(XRF.max)[3] <- "Max"
      XRF.max <- na.omit(XRF.max)                 # iniore les valeurs maquantes
      
      XRF <- merge(XRF, XRF.min, by = c("Plot.x", "variable"), all = T)
      XRF <- merge(XRF, XRF.max, by = c("Plot.x", "variable"), all = T)
    }
  }
  else{n.XRF <- 0}
  
  #### Isotope ####
  if(is.null(Isotope) == F){
    #### Selection des proxy ####
    if(is.null(Remove.Isotope) == F){Isotope <- Isotope[setdiff(names(Isotope), Remove.Isotope)]}
    if(is.null(Keep.Isotope) == F){Isotope <- Isotope[Keep.Isotope]}
    n.Iso <- ncol(Isotope)-1
    #### Mise en forme de la matrice ####
    #if(is.null(Ratio) == F){Isotope <- cbind(Isotope, Mratio)}
    Isotope <- reshape2::melt(Isotope, id = Plot.x)
    names(Isotope)[1] <- "Plot.x"
    Isotope <- na.omit(Isotope)                 # iniore les valeurs maquantes
  }
  else{n.Iso <- 0}
  
  #### MS ####
  if(is.null(MS) == F){
    MS.keep <- MS[names(MS) %in% Possible.x]
    MS <- MS[!names(MS) %in% setdiff(Possible.x, Plot.x)]
    if(is.null(Remove.MS) == F){MS <- MS[setdiff(names(MS), Remove.MS)]}
    if(is.null(Keep.MS) == F){MS <- MS[c(Keep.MS, Plot.x)]}
    n.MS <- ncol(MS)-1
    MS <- reshape2::melt(MS, id = Plot.x)
    names(MS)[1] <- "Plot.x"
    MS <- na.omit(MS)                 # iniore les valeurs maquantes
  }
  else{n.MS <- 0}
  
  #### GDGT ####
  if(is.null(GDGT) == F){
    #### Selection des proxy ####
    GDGT.keep <- GDGT[names(GDGT) %in% Possible.x]
    GDGT <- GDGT[!names(GDGT) %in% setdiff(Possible.x, Plot.x)]
    if(is.null(Remove.GDGT) == F){GDGT <- GDGT[setdiff(names(GDGT), Remove.GDGT)]}
    if(is.null(Keep.GDGT) == F){GDGT <- GDGT[c(Keep.GDGT, Plot.x)]}    
    if(is.null(Remove.GDGT) == T & is.null(Keep.GDGT) == T){
      print("Too much GDGT indexe to plot.")
      Keep.GDGT = c("BIT", "Ib.Ia", "IR6_7Me", "pCren", "CI", "MBTp5Me", "CBT5Me")
      GDGT <- GDGT[c(Keep.GDGT, Plot.x)]
    }
    n.GDGT <- ncol(GDGT)-1
    
    #### Intervalles ####
    if(is.null(GDGT.min) == F){GDGT.min <- GDGT.min[intersect(names(GDGT.min), names(GDGT))]}
    if(is.null(GDGT.max) == F){GDGT.max <- GDGT.max[intersect(names(GDGT.max), names(GDGT))]}
    
    #### Melt ####
    GDGT <- reshape2::melt(GDGT, id = Plot.x)
    names(GDGT)[1] <- "Plot.x"
    GDGT <- na.omit(GDGT)
    
    if(is.null(GDGT.min) == F){
      GDGT.min <-  reshape2::melt(GDGT.min, id = Plot.x)
      names(GDGT.min)[1] <- "Plot.x"
      names(GDGT.min)[3] <- "Min"
      GDGT.min <- na.omit(GDGT.min)
      GDGT <- full_join(GDGT, GDGT.min, by = c("Plot.x", "variable"))
    }
    if(is.null(GDGT.max) == F){
      GDGT.max <-  reshape2::melt(GDGT.max, id = Plot.x)
      names(GDGT.max)[1] <- "Plot.x"
      names(GDGT.max)[3] <- "Max"
      GDGT.max <- na.omit(GDGT.max)
      GDGT <- full_join(GDGT, GDGT.max, by = c("Plot.x", "variable"))
    }
  }
  else{n.GDGT <- 0}
  
  #### Merge datas ####
  # library(plyr) # plyr::rbind.fill
  if(is.null(MS) == F & is.null(XRF) == F & is.null(Isotope) == F & is.null(GDGT) == T){
    Label.proxy = c("XRF", "MS", "Isotopes")
    Mgeoch.core <- plyr::rbind.fill(XRF, MS, Isotope)}
  if(is.null(MS) == F & is.null(XRF) == T & is.null(Isotope) == F & is.null(GDGT) == T){
    Label.proxy = c("Isotopes", "MS")
    Mgeoch.core <- plyr::rbind.fill(Isotope, MS)}
  if(is.null(MS) == T & is.null(XRF) == F & is.null(Isotope) == F & is.null(GDGT) == T){
    Label.proxy = c("XRF", "Isotopes")
    Mgeoch.core <- plyr::rbind.fill(XRF, Isotope)}
  if(is.null(MS) == F & is.null(XRF) == F & is.null(Isotope) == T & is.null(GDGT) == T){
    Label.proxy = c("XRF", "MS")
    Mgeoch.core <- plyr::rbind.fill(XRF, MS)}
  if(is.null(MS) == F & is.null(XRF) == T & is.null(Isotope) == T & is.null(GDGT) == T){
    Label.proxy = c("MS")
    Mgeoch.core <- MS}
  if(is.null(MS) == T & is.null(XRF) == F & is.null(Isotope) == T & is.null(GDGT) == T){
    Label.proxy = c("XRF")
    Mgeoch.core <- XRF}
  if(is.null(MS) == T & is.null(XRF) == T & is.null(Isotope) == F & is.null(GDGT) == T){
    Label.proxy = c("Isotopes")
    Mgeoch.core <- Isotope}
  if(is.null(MS) == T & is.null(XRF) == T & is.null(Isotope) == T & is.null(GDGT) == F){
    Label.proxy = c("GDGT")
    Mgeoch.core <- GDGT}
  
  if(is.null(MS) == F & is.null(XRF) == F & is.null(Isotope) == F & is.null(GDGT) == F){
    Label.proxy = c("(A) XRF", "(B) MS", "(C) Isotopes and Elements Geochimistry", "(D) GDGTs")
    Mgeoch.core <- plyr::rbind.fill(XRF, MS, Isotope, GDGT)}
  if(is.null(MS) == F & is.null(XRF) == T & is.null(Isotope) == F & is.null(GDGT) == F){
    Label.proxy = c("(A) Isotopes", "(B) MS", "(C) GDGTs")
    Mgeoch.core <- plyr::rbind.fill(Isotope, MS, GDGT)}
  if(is.null(MS) == T & is.null(XRF) == F & is.null(Isotope) == F & is.null(GDGT) == F){
    Label.proxy = c("(A) XRF", "(B) Isotopes", "(C) GDGTs")
    Mgeoch.core <- plyr::rbind.fill(XRF, Isotope, GDGT)}
  if(is.null(MS) == F & is.null(XRF) == F & is.null(Isotope) == T & is.null(GDGT) == F){
    Label.proxy = c("(A) µ-XRF", "(B) MS + Spectrocolorimetry", "(C) GDGTs")
    Mgeoch.core <- plyr::rbind.fill(XRF, MS, GDGT)}
  if(is.null(MS) == F & is.null(XRF) == T & is.null(Isotope) == T & is.null(GDGT) == F){
    Label.proxy = c("(A) MS", "(B) GDGTs")
    Mgeoch.core <- plyr::rbind.fill(MS, GDGT)}
  if(is.null(MS) == T & is.null(XRF) == F & is.null(Isotope) == T & is.null(GDGT) == F){
    Label.proxy = c("(A) µ-XRF", "(B) GDGTs")
    Mgeoch.core <- plyr::rbind.fill(XRF, GDGT)}
  if(is.null(MS) == T & is.null(XRF) == T & is.null(Isotope) == F & is.null(GDGT) == F){
    Label.proxy = c("(A) Isotope", "(B) GDGTs")
    Mgeoch.core <- plyr::rbind.fill(Isotope, GDGT)}
  
  #### Categories ####
  if(is.null(Groupes) == F){
    Categorie <- as.character(Mgeoch.core$variable)
    for(i in 1:length(Groupes)){
      for(k in Groupes[[i]]){
        k <- paste("\\b", k, "\\b", sep = "")
        Categorie[grep(k, Categorie)] <- i
      }
    }
    # Categorie[str_length(Categorie) > 1] <- "0"
    Mgeoch.core <- cbind(Mgeoch.core, Categorie)
  }
  else{
    Categorie <- as.character(Mgeoch.core$variable)
    Categorie[is.factor(Mgeoch.core$variable)] <- "0"
    Mgeoch.core <- cbind(Mgeoch.core, Categorie)
    My_arrows <- NULL
  }
  
  #### Smooth curves by proxies ####
  if(is.null(Smooth.proxy) == F){
    if(Smooth.proxy == "XRF"){SubDB <- Mgeoch.core[Mgeoch.core$variable  %in% levels(XRF$variable),]}
    if(Smooth.proxy == "MS"){SubDB <- Mgeoch.core[Mgeoch.core$variable  %in% levels(MS$variable),]}
    if(Smooth.proxy == "GDGT"){SubDB <- Mgeoch.core[Mgeoch.core$variable  %in% levels(GDGT$variable),]}
    if(Smooth.proxy == "all"){SubDB <- Mgeoch.core}
    My_smooth <- geom_smooth(data = SubDB, method = "loess", se = Smooth.sd, fullrange = F, level = 0.95, linetype="solid", #formula = "x ~ y",
                             size = .8, show.legend = F, span = Smooth.param, alpha = 0.2, orientation = "y")}
  else{My_smooth <- NULL}
  
  #### CONISS  ####
  if(CONISS == T){
    Mgeoch.core.cn <- tibble::as_tibble(Mgeoch.core[,1:3])
    coniss.c <- tidypaleo::nested_data(Mgeoch.core.cn, qualifiers = Plot.x, key = variable, value = value, trans = scale, fill = 0)
    # Sous Zone
    if(is.null(Nsous.zone) == F){
      coniss.c2 <- tidypaleo::nested_chclust_coniss(coniss.c, n_groups = Nsous.zone)
      Dendro.coniss2 <- tidypaleo::layer_dendrogram(coniss.c2, aes(y = Plot.x), variable = "CONISS", sequential_facets = T, size = 0.05, alpha = 0.5)
      Line.coniss2 <- tidypaleo::layer_zone_boundaries(coniss.c2, aes(y = Plot.x), colour = "grey50", linetype = "dotted")
    }
    else{
      Dendro.coniss2 <- NULL
      coniss.c2 <- NULL
      Line.coniss2 <- NULL
    }
    
    # Zone
    if(is.null(Nzone) == F){coniss.c <- tidypaleo::nested_chclust_coniss(coniss.c, n_groups = Nzone)}
    else{coniss.c <- tidypaleo::nested_chclust_coniss(coniss.c)}
    Dendro.coniss <- tidypaleo::layer_dendrogram(coniss.c, aes(y = Plot.x), variable = "CONISS", sequential_facets = T, size = 0.05, alpha = 0.5)
    Line.coniss <- tidypaleo::layer_zone_boundaries(coniss.c, aes(y = Plot.x), colour = "grey30")
    
    if(Print.cluster.zone == T){
      df <- dplyr::select(coniss.c[[12]][[1]], Plot.x, hclust_zone)
      df <- na.omit(df)
      df <- dplyr::group_by(df, hclust_zone) 
      CONISS.bound <- dplyr::summarize(df, Min.CONISS = min(Plot.x, na.rm = T), Max.CONISS = max(Plot.x, na.rm = T))
      CONISS.bound$hclust_zone <- paste("L", rev(CONISS.bound$hclust_zone), sep = "")
      
      if(is.null(Nsous.zone) == F){
        df2 <- dplyr::select(coniss.c2[[12]][[1]], Plot.x, hclust_zone)
        df2 <- na.omit(df2)
        df2 <- dplyr::group_by(df2, hclust_zone) 
        CONISS.bound2 <- dplyr::summarize(df2, Min.CONISS = min(Plot.x, na.rm = T), Max.CONISS = max(Plot.x, na.rm = T))
        CONISS.bound2$hclust_zone <- paste("U", rev(CONISS.bound2$hclust_zone), sep = "")
        
        CONISS.bound$Category <- "Litho"
        CONISS.bound2$Category <- "Units"
        CONISS.bound.full <- rbind(CONISS.bound, CONISS.bound2)
        CONISS.bound.full <- CONISS.bound.full[order(CONISS.bound.full$Min.CONISS),]
        
        Unit <- NA
        for(i in 1:nrow(CONISS.bound.full)){
          if(grepl("L", CONISS.bound.full$hclust_zone[i]) == T){
            Unit = CONISS.bound.full$hclust_zone[i]}
          else{
            CONISS.bound.full$hclust_zone[i] <- paste("U", gsub("L", "", Unit), sep = "")
          }}
        CONISS.bound.full$hclust_zone <- ave(as.character(CONISS.bound.full$hclust_zone), CONISS.bound.full$hclust_zone, 
                                             FUN=function(x) if (length(x)>1) paste0(x[1], '.', seq_along(x)) else x[1])
        
        
        CONISS.bound <- CONISS.bound.full
      }
      
      print("CONISS layer boundaries :")
      print(CONISS.bound)
    }  
    
    if(is.null(Save.CONISS) == F){
      if(is.null(MS) == F){Fine.line <- c(match(CONISS.bound$Max.CONISS, MS.keep[[Plot.x]]), match(CONISS.bound$Min.CONISS, MS.keep[[Plot.x]]))}
      else{Fine.line <- c(match(CONISS.bound$Max.CONISS, XRF.keep[[Plot.x]]), match(CONISS.bound$Min.CONISS, XRF.keep[[Plot.x]]))}
      
      Fine.line <- sort(Fine.line[!is.na(Fine.line)])
      Fine.line <- Fine.line[1:nrow(CONISS.bound)]
      saveRDS(CONISS.bound, Save.CONISS)}
  }
  else{
    Dendro.coniss <- NULL
    Dendro.coniss2 <- NULL
    Line.coniss <- NULL
    Line.coniss2 <- NULL
    
  }
  
  #### Color settings ####
  if(is.null(Manual.color.curves)){
    A = 1:11
    Keep.col <- A[-c(6,7)] #there are 9, I exluded the two lighter hues
    if(length(unique(Mgeoch.core$Categorie)) >= 4){Keep.col <- A[-c(4,5,6,7)]}
    if(length(unique(Mgeoch.core$Categorie)) == 3){Keep.col <- A[-c(1,2,4,5,6,8,9)]}
    if(length(unique(Mgeoch.core$Categorie)) == 2){Keep.col <- A[-c(2,3,4,5,6,7,8,9,10)]}
    my_orange = brewer.pal(n = 11, "Spectral")[Keep.col] 
    orange_palette = colorRampPalette(my_orange)
    my_orange = orange_palette(length(unique(Mgeoch.core$Categorie)))
    
    if("0" %in% unique(Mgeoch.core$Categorie) == T){my_orange <- c("black", my_orange)}
    # names(my_orange) <- sort(unique(Mgeoch.core$Categorie))
  }
  else{
    my_orange <- Manual.color.curves
  }
  names(my_orange) <- unique(Mgeoch.core$Categorie)
  
  #### Graphical settings ####
  if(missing(Limites.x)){Limites.x = c(min(Mgeoch.core$Plot.x, na.rm = T), max(Mgeoch.core$Plot.x, na.rm = T))}
  Limites.zone <- c(Limites.x[1], Limites.x[2])
  Mgeoch.core <- Mgeoch.core[Mgeoch.core$Plot.x <= Limites.x[2] & Mgeoch.core$Plot.x >= Limites.x[1],]
  if(is.null(Mgeoch.core$Min) == F){Incertitude.zone <- geom_ribbon(mapping = aes(xmin = Min, xmax = Max), alpha = 0.2, size = 0.15, linetype = "dashed")}
  else{Incertitude.zone <- NULL}
  if(Display.legend == F){Display.legend <- "none"}
  if(Display.legend == T){Display.legend <- "bottom"}
  if(is.null(Axis.x.tible) == F){Title.age = Axis.x.tible}
  
  if(is.null(Nb.ticks.minors) == F){
    My_ticks <- c(Limites.zone[1], round(seq(0, Limites.zone[2], by = Nb.ticks.minors)))
    My_ticks[!grepl(paste(rep("0", str_length(Nb.ticks.minors)), collapse = ""), My_ticks)] <- ""
    My_seq <- c(Limites.zone[1], round(seq(0, Limites.zone[2], by = Nb.ticks.minors)))}
  else{
    Nb.ticks.minors = 1
    My_ticks <- c(Limites.zone[1], round(seq(0, Limites.zone[2], by = Select.interv.x)))
    My_seq <- c(Limites.zone[1], round(seq(0, Limites.zone[2], by = Select.interv.x)))}
  
  # if(max(Mgeoch.core$Plot.x, na.rm = T) > 1000){
  #   Mgeoch.core$Plot.x <- Mgeoch.core$Plot.x/10
  #   }
  if(Plot.x == "Age"){Title.age = "Time (cal. year BP)"}
  if(Plot.x == "Depth"){Title.age = "Depth (mm)"}
  if(Plot.x == "MCD"){Title.age = "Mastercore depth (mm)"}
  if(grepl("Ord", Plot.x) == T){Title.age = "Surface samples ordination"}
  
  if(is.null(Change.titles) == F){
    if(length(Label.proxy) == length(Change.titles)){Label.proxy <- Change.titles}
    else{print("You need to give the same size vector for Change.titles")}
  }
  
  if(is.null(Limites.x) == F & is.null(Manual.line) == F){Manual.line <- Manual.line[Manual.line >= Limites.x[1] & Manual.line <= Limites.x[2]]}
  if(is.null(Limites.x) == F & is.null(Manual.small.line) == F){Manual.small.line <- Manual.small.line[Manual.small.line >= Limites.x[1] & Manual.small.line <= Limites.x[2]]}
  
  #### Arrows ####
  if(Arrows == T & is.null(Groupes) == F){
    if(is.null(XRF.SD) == T & is.null(GDGT.max) == T & is.null(GDGT.min) == T){My_arrows <- unique(Mgeoch.core[c(2,4)])}
    else{My_arrows <- unique(Mgeoch.core[c(2,6)])}
    
    My_arrows <- My_arrows[match(My_arrows$variable, levels(My_arrows$variable)),]
    levels(My_arrows$variable) <- factor(My_arrows$variable)
    Max.by.kat <- aggregate(Mgeoch.core$value, by = list(Mgeoch.core$variable), FUN = max, na.rm = T)[[c(2)]]
    Min.by.kat <- aggregate(Mgeoch.core$value, by = list(Mgeoch.core$variable), FUN = min, na.rm = T)[[c(2)]]
    
    if(is.null(Arrow.y) == T){Arrow.y <- 0.015}
    My_arrows$x <- (Limites.zone[2] + Arrow.y*Limites.zone[2])
    
    My_arrows$Max <- Max.by.kat
    My_arrows$Min <- Min.by.kat
    My_arrows$Mil <- (Min.by.kat + Max.by.kat)/2
    My_arrows$Q3 <- Min.by.kat*.2 + Max.by.kat*.8
    My_arrows$Q1 <- Min.by.kat*.8 + Max.by.kat*.2
    
    # My_arrows$variable[My_arrows$Categorie == 0] <- NA
    My_arrows <- na.omit(My_arrows)
    
    My_arrows$Label.proxy <- as.character(My_arrows$variable)
    for(i in 1:nrow(My_arrows)){
      for(j in 1:length(Groupes)){
        if(My_arrows$Label.proxy[i] %in% Groupes[[j]]){
          if(is.null(names(Groupes[j])) == F){
            My_arrows$Label.proxy[i] <- names(Groupes[j])}}}}
    
    # My_arrows$Label.proxy <- gsub("\\.", " ", My_arrows$Label.proxy)
    My_arrows$Label.proxy <- gsub("\\.", "\n", My_arrows$Label.proxy)
    if(is.null(Reverse.arrow) == F){
      for(i in Reverse.arrow){
        New_Q1 <- My_arrows[i,"Q3"]
        New_Q3 <- My_arrows[i,"Q1"]
        My_arrows[i,"Q1"] <- New_Q1
        My_arrows[i,"Q3"] <- New_Q3
      }
    }
    My_arrows$Arrow.y <- Arrow.y
    My_Plab <- geom_text(data = My_arrows, inherit.aes = F,
                         aes(x = Mil,
                             y = (x + 1.5*Arrow.y*x),
                             label = Label.proxy#,
                             # color = Categorie
                         )
    )
    
    My_arrows.draw <- geom_segment(data = My_arrows, inherit.aes = F,
                                   aes(x = Q1,
                                       y = x,
                                       xend = Q3,
                                       yend = x, color = Categorie),
                                   size = .6, arrow = arrow(length = unit(0.3, "cm")))
    
  }
  else{
    My_Plab <- NULL
    My_arrows.draw <- NULL
  }
  
  #### Climate zone names ####
  if(is.null(Name.zone) == F & is.null(Zone.clim) == F){
    library("tableHTML")
    library("ggh4x")
    if(is.null(Keep.zone.clim.full) == F){Zone.clim <- Keep.zone.clim.full}
    
    data_hline <- data.frame(variable = as.factor(c(levels(Mgeoch.core$variable), rep(" ", length(Name.zone)))),
                             SECT.NUM = c(rep(NA,nlevels(Mgeoch.core$variable)), Name.zone),
                             MY_col = c(rep(NA,nlevels(Mgeoch.core$variable)), Temp.zone),
                             Lab.size = c(rep(NA, nlevels(Mgeoch.core$variable)), rep(5, length(Name.zone))),
                             xmin = c(rep(NA, nlevels(Mgeoch.core$variable)), Zone.clim[odd(seq(1, length(Zone.clim)))]),
                             xmax = c(rep(NA, nlevels(Mgeoch.core$variable)), Zone.clim[even(seq(1, length(Zone.clim)))]),
                             ymin = c(rep(NA, nlevels(Mgeoch.core$variable)), rep(4, length(Name.zone))),
                             ymax = c(rep(NA, nlevels(Mgeoch.core$variable)), rep(8, length(Name.zone))))
    
    if(is.null(Name.zone2) == F & is.null(Zone.clim2) == F){
      data_hline2 <- data.frame(variable = as.factor(rep(" ", length(Name.zone2))), SECT.NUM = Name.zone2, MY_col = Temp.zone2,
                                xmin = Zone.clim2[odd(seq(1, length(Zone.clim2)))],
                                xmax = Zone.clim2[even(seq(1, length(Zone.clim2)))], Lab.size = rep(4.95, length(Name.zone2)),
                                ymin = rep(0, length(Name.zone2)), ymax = rep(4, length(Name.zone2)))
      
      data_hline <- rbind(data_hline, data_hline2)
      if(is.null(Panel.param) == T){Panel.param = 0.5}
    }
    else{if(is.null(Panel.param) == T){Panel.param = 0.35}}
    
    data_hline <- data_hline[data_hline$SECT.NUM %in% data_hline$SECT.NUM[data_hline$xmax >= Limites.x[1] & data_hline$xmin <= Limites.x[2]],]
    data_hline$xmin[which.min(data_hline$xmin)] <- Limites.x[1]
    data_hline$xmax[which.max(data_hline$xmax)] <- Limites.x[2]
    
    Panel.size <- ggh4x::force_panelsizes(cols = c(rep(1, nlevels(Mgeoch.core$variable)), Panel.param))
    Rect_core <- geom_rect(data = data_hline, inherit.aes = F,
                           mapping = aes(xmin = ymin, xmax = ymax, ymin = xmin, ymax = xmax, fill = MY_col, color = MY_col),
                           alpha=0.6, #color = "grey30",
                           # linewidth = .5, show.legend = F, 
                           na.rm = T
    )
    Text_core <- geom_text(data = data_hline, inherit.aes = F,
                           aes(x = (ymin+ymax)/2, y = (xmin+xmax)/2, label = SECT.NUM, size = Lab.size),
                           vjust = 0.5, hjust = 0.5, color = "grey20" , na.rm = T)}
  else{
    Panel.size <- NULL
    Rect_core <- NULL
    Text_core <- NULL
  }
  
  #### C14 dates ####
  if(is.null(C14.path) == F){
    C14.table <- data.frame(read.csv(C14.path))
    if(is.null(Limites.x) == F ){C14.table <- C14.table[C14.table$depth >= Limites.x[1]/10 & C14.table$depth <= Limites.x[2]/10,]}
    
    library(Bchron)
    for(i in 1:nrow(C14.table)){
      if(C14.table$age[i] > 0){
        singleResult <- BchronCalibrate(C14.table$age[i], C14.table$error[i], "intcal20", allowOutside = T)
        
        C14.table$age[i] <- round(mean(singleResult$Date1$ageGrid), digits = 0)
      }
    }
    C14.table$age[C14.table$age > 1000] <- paste(round(C14.table$age[C14.table$age > 1000], digits = -2)/1000, "k")
    if(any(names(C14.table) == "labID") == T){
      data_C14 <- data.frame(variable = as.factor(c(levels(Mgeoch.core$variable), rep("C14", 2*nrow(C14.table)))),
                             SECT.NUM = c(rep(NA,nlevels(Mgeoch.core$variable)), C14.table$labID, C14.table$labID),
                             Lab.etiquette = c(rep(NA, nlevels(Mgeoch.core$variable)), rep(NA, nrow(C14.table)), C14.table$age),
                             x = c(rep(NA, nlevels(Mgeoch.core$variable)), rep(NA, nrow(C14.table)), C14.table$depth*10),
                             y = c(rep(NA, nlevels(Mgeoch.core$variable)), rep(8, length(C14.table$depth)), rep(0, length(C14.table$depth)))
      )
      
      if(is.null(Panel.param) == T){Panel.param = 0.35}
      Panel.size2 <- ggh4x::force_panelsizes(cols = c(rep(1, nlevels(Mgeoch.core$variable)), Panel.param))
      
      pC14 <- geom_point(data = data_C14, inherit.aes = F,
                         aes(x = y, y = x), color = "black",
                         size = 2, alpha = 1, na.rm = T)
      
      pC14_txt <- geom_text_repel(data = data_C14, inherit.aes = F,
                                  mapping = aes(x = y, y = x, label = Lab.etiquette),
                                  nudge_y = 10, segment.inflect = F, segment.square = F, segment.curvature = 0.2,
                                  force = 5, nudge_x  = 200, direction = "y", hjust = 0, na.rm = T,
                                  size = 4, parse = F, segment.size = 0.5, segment.colour = "grey50"
      )
    }
    else{
      print("The table of C14 must have its first column as labID.")
      Panel.size2 <- NULL; pC14 <- NULL; pC14_txt <- NULL
    }
  }
  else{
    Panel.size2 <- NULL; pC14 <- NULL; pC14_txt <- NULL
  }
  
  #### Plot ####
  p1 <- ggplot(data = Mgeoch.core, mapping =  aes(y = Plot.x,
                                                  x = value,
                                                  color = Categorie, fill = Categorie)) + 
    #### Geoch ####
  tidypaleo::geom_lineh(linewidth = 0.1) +
    geom_hline(yintercept = rep(Manual.line, nlevels(Mgeoch.core$variable)), col = "grey10", lty = 2, alpha = 0.6)+
    geom_hline(yintercept = rep(Manual.small.line, nlevels(Mgeoch.core$variable)), col = "grey30", linetype = "dotted", alpha = 0.6)+
    Incertitude.zone +
    geom_point(size = 0.8, na.rm = T) +
    My_smooth +
    My_arrows.draw +
    My_Plab + 
    scale_y_reverse(breaks = My_seq, labels = My_ticks, expand = c(0.05, 0)) +
    scale_x_continuous(name = NULL, position = "top",  breaks = scales::pretty_breaks(n = 4), limits = NULL) +
    # scale_x_continuous(name = NULL, position = "top",  breaks = scales::pretty_breaks(n = 4), limits = NULL, expand = c(0,0.1)) +
    labs(x = NULL, y = Title.age) +
    Hiatus.area +
    tidypaleo::facet_geochem_gridh(vars(variable), renamers = New.lab, units = My.unit, default_units = "%") +
    Panel.size + Dendro.coniss + Line.coniss + Line.coniss2 +
    scale_fill_manual(values = my_orange, guide = "none")+
    
    #### Zone clim ####
  new_scale_fill()+
    Zone.clim.to.plot +
    scale_fill_manual(values = rev(Rect.color.scale), guide = "none", name = NULL, labels = NULL, breaks = NULL, na.translate = FALSE)+
    scale_size_continuous(range = c(3.5,5))+
    Rect_core + Text_core +
    
    guides(color = guide_legend(nrow = 1))+        # force les legendes a salligner sur une unique ligne
    scale_color_manual(values = my_orange, name = "Proxies: ", labels = names(Groupes))+
    Panel.size2 + pC14_txt + pC14 +
    #### Theme ####
  theme(
    axis.title.x=element_text(size=13),              # XRF
    axis.text.x.top = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    axis.title.y=element_text(size=15),              # AGE
    axis.line.x = element_line(lineend = "butt", ),                   # fait apparaitre seulement l'axe des y
    axis.line.y = element_line(lineend = "butt"),    # fait apparaitre seulement l'axe des x, coupe au bout
    legend.position = Display.legend,                        # permet de mettre le carre des legendes en bas
    legend.key = element_blank(),  # carré autours du symbole
    legend.direction = "horizontal", plot.background = element_blank(), panel.grid = element_blank(),
    panel.spacing = unit(0.04, "lines"),
    panel.background=element_blank(), strip.background = element_blank(),
    strip.placement = "outside", panel.border = element_blank(), strip.clip = "off",
    strip.text = element_text(size = 11),
    legend.justification = c("left"),               # left, top, right, bottom
    plot.margin = unit(c(0,0,0,0), "lines")
  )
  
  #### Grob ####
  if(length(which(c(is.null(XRF), is.null(MS), is.null(GDGT), is.null(Isotope)) == F)) > 1){
    z <- ggplot2::ggplotGrob(p1)
    my_pos = 8 # Il faut essyer au hasard de changer cette valeur si pb apres maj
    z <- gtable::gtable_add_rows(z, z$height[my_pos], pos = c(my_pos-1))  # New row added below row my_pos (avant 6)
    # z <- gtable::gtable_add_rows(z, unit(10/10, "cm"), 14)
    
    N.tot <- c(n.XRF, n.XRF+n.MS, n.XRF+n.MS+n.Iso, n.XRF+n.MS+n.Iso+n.GDGT)
    N.tot <- unique(N.tot)
    N.tot <- (N.tot*2)
    N.tot <- c(1, N.tot)
    N.tot <- N.tot + 5
    for(i in 1:length(Label.proxy)){
      z <- gtable::gtable_add_grob(z,
                                   list(grid::rectGrob(gp = grid::gpar(col = NA, fill = NA, size = .5)),
                                        grid::textGrob(Label.proxy[i], gp = grid::gpar(cex = 1.2, fontface = 'bold', col = "black"))),
                                   t = my_pos, l=N.tot[i], b = my_pos, r=N.tot[i+1], name = c("a", "b"))
      z <- gtable::gtable_add_cols(z, unit(2/10, "line"), N.tot[i+1])
      N.tot <- N.tot + 1
    }
    Ptot <- z
  }
  else{Ptot <- p1}
  
  #### Save plot and export ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(Ptot, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.plot)}}
  else{print(Ptot)}
  
  if(Print.cluster.zone == T & CONISS == T){return(CONISS.bound)}
  if(Print.cluster.zone == F & CONISS == F){
    if(return.plot == T){
      return(Ptot);print("yes")
    }
    else{return(Mgeoch.core)}
  }
}

PCA.XRF <- function(MP, transp_OK, Site.name, Type.samples, Ellipse, Show.centroid, Show.color, Show.lab.PCA = T,
                    Csv.sep, Scale.PCA, Groupes, Cluster.core, Cluster.core.lab, return.plot, 
                    Legends.pos = "right", Leg.nrows = NULL, Color.vectors = NULL, Alpha.ellipse = 0.2,
                    X.pos, Y.pos, GDGT, Color.choice, Save.path, Manu.lim, Save.plot, H, W){
  #### Settings ####
  library(vegan)
  library(FactoMineR)
  library(factoextra)
  if(missing(Show.color)){Show.color = F}
  if(missing(return.plot)){return.plot = F}
  if(missing(Show.centroid)){Show.centroid = F}
  if(Show.centroid == F){Centroide = "quali"}
  if(Show.centroid == T){Centroide = NULL}
  if(missing(Csv.sep)){Csv.sep = "\t"}
  if(missing(Cluster.core)){Cluster.core = "Age"}
  if(missing(Cluster.core.lab)){Cluster.core.lab = "Age"}
  if(missing(Scale.PCA)){Scale.PCA = 1}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Ellipse)){Ellipse = F}
  if(missing(transp_OK)){transp_OK = T}
  if(missing(Groupes)){Groupes = NULL}
  if(missing(Manu.lim)){Manu.lim = NULL}
  if(missing(Type.samples)){Type.samples = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(GDGT)){GDGT = F}
  if(missing(Site.name)){Site.name = NULL}
  if(missing(Color.choice)){Color.choice = NULL}
  if(missing(X.pos)){X.pos = "bottom"}
  if(missing(Y.pos)){Y.pos = "left"}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Data pulishing ####
  Remove.name <- c("Top", "Bottom", "Depth", "Age", "Litho", "Units", "MCD")
  Keep.xdata <- MP[intersect(names(MP), Remove.name)]
  MP <- MP[setdiff(names(MP), Remove.name)]
  #### Pour les GDGTs ####
  if(GDGT == T){
    names(MP) <- gsub("f.", "", names(MP))
    names(MP) <- gsub("_5Me", "", names(MP))
    names(MP) <- gsub("_6Me", "\\'", names(MP))
    names(MP) <- gsub("_7Me", "\\''", names(MP))
  }
  
  #### Groupes ####
  if(is.null(Groupes) == F){
    Groupes <- reshape2::melt(Groupes)
    Groupes <- Groupes$L1[match(names(MP), Groupes$value)]
    Groupes[is.na(Groupes)] <- "Unknown"
    Groupes <- factor(Groupes)
  }
  else{Groupes <- "darkgrey"}
  
  if(is.null(Leg.nrows) == T){
    if(nlevels(Groupes) > 8){Nb.columns = 2}
    else{Nb.columns = 1}
  }
  else{
    Nb.columns <- NULL
  }
  
  #### Fullfill NAs ####
  if(nlevels(as.factor(is.na(MP))) >= 2){
    library(missMDA)
    nb <- estim_ncpPCA(MP, ncp.max=5) ## Time consuming, nb = 2
    MP.comp <- imputePCA(MP, ncp = nb[[1]])
    MP <- MP.comp$completeObs
  }
  
  #### Transforming the data + PCA calcul ####
  if(transp_OK == FALSE){MP.pca <- rda(MP, scale = T)}
  else{MP.pca <- rda(log1p(MP), scale = T)}
  site.score  <- scores(MP.pca, choices=c(1,2), display = "wa", scaling = 2)         # valeur PCA age
  taxa.score <- scores(MP.pca, choices=c(1,2), display = "species", scaling = 2)
  
  PC1.MP <- round(MP.pca[["CA"]][["eig"]][["PC1"]]/MP.pca[["tot.chi"]]*100, digits = 0)              # calcul PC1
  PC2.MP <- round(MP.pca[["CA"]][["eig"]][["PC2"]]/MP.pca[["tot.chi"]]*100, digits = 0)              # calcul PC2
  MP.pca <- PCA(MP, graph = FALSE, scale.unit = transp_OK)
  
  #### Labels ####
  A = round(MP.pca$eig[1,2], digits = 0)
  B = round(MP.pca$eig[2,2], digits = 0)
  
  Xlab <- labs(x = substitute(paste("PCA"[1], ~ "(", A, " %)", sep = " " )),
               y = substitute(paste("PCA"[2], ~ "(", B, " %)", sep = " " )))
  Axes = c(1,2)
  
  if(is.null(Site.name) == F & is.null(Type.samples) == F){My_title <- paste("PCA", Site.name, Type.samples, sep = " ")}
  if(is.null(Site.name) == T & is.null(Type.samples) == F){My_title <- paste("PCA", Site.name, sep = " ")}
  if(is.null(Site.name) == F & is.null(Type.samples) == T){My_title <- paste("PCA", Type.samples, sep = " ")}
  if(is.null(Site.name) == T & is.null(Type.samples) == T){My_title <- NULL}
  if(Show.lab.PCA == F){My_title <- gsub("PCA ", "", My_title)}
  
  #### Color settings ####
  if(is.null(Color.vectors) == T){my_orange <- c("#A6CEE3", "#EB4445", "royalblue", "#87CA6A", "darkred", "#F8A361", "grey", "pink", "grey10", "darkorange")}
  else{my_orange <- Color.vectors}
  
  names(my_orange) <- sort(unique(Groupes))
  Keep.col2 <- 1:11
  if(is.numeric(Keep.xdata[[Cluster.core]]) == T){
    my_orange2 = brewer.pal(n = 11, "RdYlBu")[Keep.col2[-c(3,4,5,7,8,9)]] 
    orange_palette2 = colorRampPalette(my_orange2)
    my_orange2 = rev(orange_palette2(length(seq(min(Keep.xdata[[Cluster.core]]), max(Keep.xdata[[Cluster.core]]), by = 200))))
    Scale.fill <- scale_fill_gradientn(colours = my_orange2, guide = "colourbar", 
                                       name = Cluster.core.lab,
                                       breaks = seq(round(min(Keep.xdata[[Cluster.core]]),digits = -3), max(Keep.xdata[[Cluster.core]]), by = 1000),
                                       na.value = "white")
    if(Ellipse == T){
      print("Be carefull ! It is not possible to have ellipse with numerical core clustering.")
      Ellipse = F}
  }
  else{
    
    if(is.null(Color.choice) == F){my_orange2 <- rev(Color.choice)}
    else{my_orange2 = brewer.pal(n = 11, "RdYlBu")[Keep.col2[-c(2,4,5,7,8,10)]]}
    
    orange_palette2 = colorRampPalette(my_orange2)
    my_orange2 = rev(orange_palette2(length(unique(Keep.xdata[[Cluster.core]]))))
    if(Show.color == T){
      print("My color scale is :")
      print(my_orange2)
    }
    
    Scale.fill <- scale_fill_manual(values = my_orange2, name = Cluster.core.lab)
  }
  
  #### PLOT ####
  p <- fviz_pca_biplot(MP.pca, axes = Axes,
                       geom.ind = "point",
                       pointshape = 21,
                       pointsize = 2.1,
                       fill.ind = Keep.xdata[[Cluster.core]],
                       alpha.ind = 0.8,
                       col.ind = "black",
                       repel = T,
                       invisible = Centroide, # enlève ou ajoute le centroïde
                       addEllipses = Ellipse, ellipse.level = 0.9, ellipse.type = "convex", # convex or norm
                       ellipse.alpha = Alpha.ellipse,
                       title = My_title,
                       col.var = Groupes) +
    Xlab+ Scale.fill + 
    guides(colour = guide_legend(override.aes=list(size=0.5), nrow = Leg.nrows), fill = guide_legend(ncol = Nb.columns, nrow = Leg.nrows))+
    scale_x_continuous(position = X.pos) + 
    scale_y_continuous(position = Y.pos) +
    scale_color_manual(values = my_orange, name = "Proxies")+
    theme(plot.background = element_blank(), legend.position = Legends.pos,
          panel.border = element_rect(NA, "black", linewidth = 1),
          plot.margin=unit(c(0,0,0,0),"cm"),
          panel.grid = element_line(linetype = "dashed"),
          axis.line = element_blank())
  
  #### Export data ####
  PCA <- data.frame(MP.pca$ind$coord)
  PCA <- PCA[,1:2]
  names(PCA) <- c("PC1", "PC2")
  
  if(is.null(Save.path) == F){
    Site.name <- gsub(" ","_",Site.name)
    Save.path.Site <- gsub("\\.csv", "_PCA_core.csv", Save.path)
    Save.path.Taxon <- gsub("\\.csv", "_PCA_elemt.csv", Save.path)
    write.table(PCA, file = Save.path.Site, col.names = T, sep = ",", row.names = F)
    write.table(t(taxa.score), file = Save.path.Taxon, col.names = TRUE, sep = ",")
  }
  if(is.null(Save.plot) == F){
    print(p)
    dev.off()}
  if(return.plot == T){return(p)}
  else{return(PCA)}
}

Plot.BIT.IIIa <- function(MG, Select.param, Cluster.core, Manu.lim, Vline, Hline, Leg.pos, Remove.Cluster.from.reg.lin = NULL,
                          No.lines = F, Add.reg.lin = F, R2.pos = NULL, Symbol.path = NULL, Symbol.pos = NULL, 
                          Bisectrice = F, Ellipse = F, Alpha.ellipse = 0.2, RL.formula = "Linear",
                          Title, Color.choice, Cluster.core.lab, W, H, Save.plot){
  #### Settings ####
  if(missing(H)){H = NULL}
  if(missing(W)){W = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Select.param)){Select.param = NULL}
  if(missing(Cluster.core)){Cluster.core = NULL}
  if(missing(Cluster.core.lab)){Cluster.core.lab = "Lithology"}
  if(missing(Manu.lim)){Manu.lim = NULL}
  if(missing(Vline)){Vline = NULL}
  if(missing(Leg.pos)){Leg.pos = NULL}
  if(missing(Hline)){Hline = NULL}
  if(missing(Color.choice)){Color.choice = NULL}
  if(missing(Title)){Title = NULL}
  
  if(is.null(Select.param) == T){Check.origin <- MG[,c("BIT","IIIa.IIa")]}
  else{Check.origin <- MG[,Select.param]}
  Lab.y <- names(Check.origin)[1]
  Lab.x <- names(Check.origin)[2]
  names(Check.origin) <- c("y", "x")
  
  #### Cluster core ####
  if(is.null(Cluster.core) == F){
    Check.origin <- cbind(Check.origin, MG[Cluster.core])
    names(Check.origin)[3] <- "cluster"
    
    #### Color settings ####
    my_orange <- c("#A6CEE3", "#EB4445", "royalblue", "#87CA6A", "darkred", "#F8A361", "grey")
    Keep.col2 <- 1:11
    if(is.numeric(Check.origin[[Cluster.core]]) == T){
      my_orange2 = brewer.pal(n = 11, "RdYlBu")[Keep.col2[-c(3,4,5,7,8,9)]]
      orange_palette2 = colorRampPalette(my_orange2)
      my_orange2 = rev(orange_palette2(length(seq(min(Check.origin[[Cluster.core]]), max(Check.origin[[Cluster.core]]), by = 200))))
      Scale.fill <- scale_color_gradientn(colours = my_orange2, guide = "colourbar",
                                          name = Cluster.core.lab,
                                          breaks = seq(round(min(Check.origin[[Cluster.core]]),digits = -3), max(Check.origin[[Cluster.core]]), by = 1000),
                                          na.value = "white")}
    else{
      if(is.null(Color.choice) == F){my_orange2 <- rev(Color.choice)}
      else{my_orange2 = brewer.pal(n = 11, "RdYlBu")[Keep.col2[-c(2,4,5,7,8,10)]]}
      orange_palette2 = colorRampPalette(my_orange2)
      my_orange2 = rev(orange_palette2(length(unique(Check.origin$cluster))))
      Scale.fill <- scale_color_manual(values = my_orange2, name = Cluster.core.lab)
    }
    
  }
  
  #### Param ####
  if(is.null(Leg.pos) == F){Title <- ggtitle(Title)}
  if(is.null(Leg.pos) == T){Leg.pos <- "left"}
  if(is.null(Manu.lim) == F){My_lims <- lims(x = c(Manu.lim[1],Manu.lim[2]),  y = c(Manu.lim[3],Manu.lim[4]))}
  else{My_lims <- NULL}
  
  if(is.null(Vline) == T){
    V1 <- geom_vline(xintercept = 0.86,lty = "dotted", col = "darkorange", linewidth= .8)
    V2 <- geom_vline(xintercept = 1.2,lty = "dotted", col = "royalblue", linewidth= .8)
  }
  else{
    V1 <- geom_vline(xintercept = Vline[1],lty = "dotted", col = "darkorange", linewidth= .8)
    V2 <- geom_vline(xintercept = Vline[2],lty = "dotted", col = "royalblue", linewidth= .8)
  }
  
  if(is.null(Hline) == T){
    H1 <- NULL
    H2 <- geom_hline(yintercept = 0.5,lty = "dotted", col = "royalblue", linewidth= .8)}
  else{
    H1 <- geom_hline(yintercept = Hline[1], lty = "dotted", col = "darkorange", linewidth= .8)
    H2 <- geom_hline(yintercept = Hline[2], lty = "dotted", col = "royalblue", linewidth= .8)}
  
  if(No.lines == T){V1 = NULL; V2 = NULL; H1 = NULL; H2 = NULL}
  # Check.origin$Lab.etiquette[Mtot.melt$Lab.etiquette == "Rib.Index2"] <- "Index[2]"
  
  if(Lab.x == "IIIa.IIa"){Lab.x <- "IIIa/IIa"}
  if(Lab.y == "Ib.Ia"){Lab.y <- "Ib/Ia"}
  if(Lab.y == "BIT"){Lab.y <- "BIT index"}
  if(Lab.x == "pCren"){Lab.x <- "%(Cren)"}
  if(Lab.y == "MBTp5Me"){Lab.y <- expression(MBT*minute[5~Me])}
  if(Lab.x == "MBTp6Me"){Lab.x <- expression(MBT*minute[6~Me])}
  
  if(Add.reg.lin == T){
    if(RL.formula == "Linear"){my_lm <- geom_smooth(method = "lm", se = F, colour = "grey30", linewidth = .5)}
    if(RL.formula == "Hyperbolic"){
      if(is.null(Remove.Cluster.from.reg.lin) == F){DF.RL <- Check.origin[!Check.origin$cluster %in% Remove.Cluster.from.reg.lin,]}
      else{DF.RL <- Check.origin}
      
      # fit <- nls(pCren ~ a / (MI + b) + c, data = GDGT.Zal, start = list(a = 1, b = 0.01, c = 0.1))
      fit <- minpack.lm::nlsLM(pCren ~ a / (MI + b) + c, data = GDGT.Zal, start = list(a = 1, b = 0.01, c = 0.1))
      pred <- predict(fit)
      RSS <- sum((GDGT.Zal$pCren - pred)^2)
      TSS <- sum((GDGT.Zal$pCren - mean(GDGT.Zal$pCren))^2)
      R2 <- 1 - RSS/TSS
      R2_text <- paste0("R² = ", round(R2, 2))
      
      summary_stats <- summary(fit)
      coefs <- summary_stats$coefficients
      t_vals <- coefs[, "Estimate"] / coefs[, "Std. Error"]
      p_vals <- 2 * pt(abs(t_vals), df = summary_stats$df[2], lower.tail = FALSE)
      if(as.numeric(p_vals[1]) < 0.001){R2_text <- paste0("R² = ", round(R2, 2), ", p > 0.001")}
      
      hyperbolic_model <- function(formula, data, ...){
        # nls(formula, data = data, start = list(a = 1, b = 0.01, c = 0.1)))
        minpack.lm::nlsLM(formula, data = data, start = list(a = 1, b = 0.01, c = 0.1))
      }
      my_lm <- stat_smooth(
        data = DF.RL,
        method = hyperbolic_model, formula = y ~ a / (x + b) + c,
        method.args = list(start = list(a = 1, b = 0.01, c = 0.1)),
        se = FALSE, colour = "grey30", linewidth = .5)
    }
  }    
  else{my_lm <- NULL}
  
  if(Bisectrice == T){my_bis <- geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70")}
  else{my_bis <- NULL}
  
  #### Add R2 ####
  if(is.null(R2.pos) == F){
    if(RL.formula == "Hyperbolic"){
      if(R2.pos == "bottomleft"){R2.y = "bottom"; R2.x = "left"}
      if(R2.pos == "bottomright"){R2.y = "bottom"; R2.x = "right"}
      if(R2.pos == "topright"){R2.y = Inf; R2.x = Inf}
      if(R2.pos == "topleft"){R2.y = "top"; R2.x = "left"}
      if(R2.pos == "none"){R2.y = "none"; R2.x = "none"}
      
      Add.r2 <- annotate("text", x = R2.x, y = R2.y, label = R2_text, hjust = 1.1, vjust = 2, size = 3)
    }
    if(RL.formula == "Linear"){
      if(R2.pos == "bottomleft"){R2.y = "bottom"; R2.x = "left"}
      if(R2.pos == "bottomright"){R2.y = "bottom"; R2.x = "right"}
      if(R2.pos == "topright"){R2.y = "top"; R2.x = "right"}
      if(R2.pos == "topleft"){R2.y = "top"; R2.x = "left"}
      if(R2.pos == "none"){R2.y = "none"; R2.x = "none"}
      
      Add.r2 <- stat_poly_eq(label.y = R2.y, label.x = R2.x, 
                             size = 3, small.r = F, vstep = 0.07, p.digits = 3, na.rm = T,  
                             aes(label =  sprintf("%s*\", \"*%s" ,
                                                  after_stat(rr.label),
                                                  # after_stat(r.squared),
                                                  after_stat(p.value.label)
                             )))}
  }
  else{Add.r2 <- NULL}
  
  #### Ajout symbole ####
  if(is.null(Symbol.path) == F){
    if(is.null(Symbol.pos)== T){Symbol.pos <- c(.9, .9, .16)}
    if(grepl("\\.png", Symbol.path)){
      library(png)
      library(grid)
      img <- readPNG(Symbol.path)
      g <- rasterGrob(x = Symbol.pos[1], y = Symbol.pos[2], width = Symbol.pos[3], height = Symbol.pos[3], img, interpolate = T)
    }
    
    if(grepl("\\.xml", Symbol.path)){
      library(grImport)
      img <- readPicture(Symbol.path)
      g <- pictureGrob(x = Symbol.pos[1], y = Symbol.pos[2], width = Symbol.pos[3], height = Symbol.pos[3], img)
    }
    
    Logo <- annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  }
  else{Logo <- NULL}
  
  #### Ajout ellipse / area ####
  if(Ellipse == T){
    library(tidyverse)
    Data_area <- Check.origin %>%
      group_by(cluster) %>%
      slice(chull(x, y))
    
    my_Area <- geom_polygon(data = Data_area, mapping = aes(fill = cluster), alpha = Alpha.ellipse)
    Scale.fill.area <- scale_fill_manual(values = my_orange2, name = Cluster.core.lab)
    
  }
  else{my_Area <- NULL; Scale.fill.area <- NULL}
  
  #### ggplot ####
  p <- ggplot(Check.origin, aes(y = y, x = x))+
    my_bis + my_Area + Scale.fill + Scale.fill.area + My_lims + labs(x= Lab.x, y = Lab.y)+
    geom_point(mapping = aes(color = cluster), shape = 16, size = 2, alpha = 0.8, na.rm = T)+
    V1 + V2 + H1 + H2 + Title + my_lm + Add.r2 + Logo + 
    theme_bw()+theme(legend.position = Leg.pos, plot.background = element_blank(),
                     panel.border = element_rect(NA, "black", linewidth = 1),
                     panel.grid = element_line(linetype = "dashed"),
                     plot.margin=unit(c(0,0,0,0.1),"cm"), panel.background = element_blank(), legend.background = element_blank(),
                     axis.line = element_blank())
  
  print(p)
  
  #### Save plot and export ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(p, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.plot)}}
  return(p)
}

DP_plot_algue <- function(MP_lake, MA_lake, Plot.x = "Age", Save.plot, H, W, Time.window = 200,
                          nlake, CONISS, Nzone){
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(nlake)){nlake = NULL}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  ### Prepa data
  Prof <- MA_lake[1,]
  Prof <- Prof[colnames(MP_lake)]
  Age <- MA_lake[Plot.x,]
  Age <- Age[colnames(MP_lake)]
  
  ### Fusion taxons rares
  Max_seuil = 0.01
  MP_max <- MP_lake[0]
  MP_max[, "max"] <- apply(MP_lake[,], 1, max)
  MP_dom <- MP_lake[rowSums(MP_max)>Max_seuil,]
  MP_rare <- MP_lake[rowSums(MP_max)<=Max_seuil,]
  MP_dom["Other algae" ,] <- colSums(MP_rare)
  
  ## Mongolie
  # InfoPol  <- read.csv(file="Import/POLLEN/Index_pollen_Mongolia.csv",sep=",",dec=".", header = T, stringsAsFactors = F, row.names = 1)
  # InfoPol  <- cbind(InfoPol, Couleur=rep(InfoPol$AP.NAP))
  # InfoPol$Couleur<- gsub('NAP', '#f27041ff', InfoPol$Couleur)
  # InfoPol$Couleur<- gsub('AP', '#6db38fff', InfoPol$Couleur)
  # InfoPol$Couleur<- gsub('NaN', '#cfce90ff', InfoPol$Couleur)
  # 
  ### Couleur
  # C <- subset(InfoPol, select=c(Couleur))
  # row.names(C)<- InfoPol$Nom
  # MP_color <- cbind(MP_dom[0], C[rownames(MP_dom),])
  # couleur = as.character(MP_color[[1]])
  
  ### Calcul AP et NAP
  #AP <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$AP.NAP == "AP"]),])
  #NAP <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$AP.NAP == "NAP"]),])
  #MP_dom <- rbind(MP_dom, AP=AP, NAP=NAP)
  
  spec <- data.frame(t(MP_dom))
  spec <- spec*100                   # On passe en pourcentage
  
  
  spec <- spec[order(as.double(Age)),]
  Age <- Age[order(as.double(Age))]
  
  #### Lake name ####
  if(is.null(nlake) == F){
    title(main = paste("Pollen diagram of", nlake, "Lake"), line = 7, cex.main = 2.7, font.main = 1, adj = 0)
    My_lake_name <- paste("Pollen diagram of", nlake, "Lake")
  }
  else{My_lake_name <- NULL}
  
  if(Plot.x == "Age"){Title.x <- "Time (yr cal BP)"}
  else{Title.x <- "Core depth (cm)"}
  
  ### CONISS Clustering
  if(CONISS == TRUE){
    diss <- dist(sqrt(t(MP_dom)/100)^2)
    Clust_DP <- chclust(diss, method = "coniss")
    #bstick(Clust_DP, 10)
  }
  else{Clust_DP = NULL}
  
  ### Plot
  p <- strat.plot(spec, yvar = as.double(Age), 
                  y.rev = TRUE, 
                  scale.percent = TRUE,     # True pour des pourcentages
                  srt.xlabel = 45,           # Rotation de 45 des noms de taxon
                  title = My_lake_name,
                  ylabel = Title.x,
                  y.tks = seq(-50,round(max(Age),digits=-2), Time.window),
                  plot.poly = TRUE, 
                  col.poly =      "#a2dbe1ff",         # bleu clair
                  col.bar =       "#8c92a1ff",         # Fait apparaitre les traits continus
                  col.poly.line = "#8c92a1ff",         # gris fonce
                  exag=FALSE,                # Fait apparaitre la zone x10
                  col.exag="auto",           # Zone x10 couleur auto
                  clust = Clust_DP,
                  clust.width=0.075 
                  #fun1=sm.fun                # application ou non du lissage
  )
  if(CONISS == TRUE & Nzone > 0){addClustZone(p, Clust_DP, Nzone, 
                                              lwd=1.5, lty=2, col="grey25") 
  }
  if(is.null(Save.plot) == F){dev.off()}
}

DP_plot_Fungal <- function(MP_lake, MA_lake, Conc.Plot, Cortege, Pol.Sum, Path.index, Save.plot, H, W, Log.trans,
                           nlake, CONISS, Nzone, Auto.lab.NPP, Time.window, Keep.cortege = NULL, Order.fungal = NULL,
                           Seuil.NPP, Diversity.pollen.ratio, Col.cortege, NPP.conc.round, Limites = NULL,
                           Spore.diversity, Richesse.type, Spore.influx, Sort.taxon, Save.path){
  #### Default values ####
  if(missing(Log.trans)){Log.trans = FALSE}
  if(missing(Col.cortege)){Col.cortege = FALSE}
  if(missing(Cortege)){Cortege = FALSE}
  if(missing(CONISS)){CONISS = FALSE}
  if(missing(Nzone)){Nzone = FALSE}
  if(missing(Conc.Plot)){Conc.Plot = FALSE}
  if(missing(Pol.Sum)){Pol.Sum = NULL}
  if(missing(Spore.diversity)){Spore.diversity = FALSE}
  if(missing(Richesse.type)){Richesse.type = FALSE}
  if(missing(Diversity.pollen.ratio)){Diversity.pollen.ratio = FALSE}
  if(missing(Auto.lab.NPP)){Auto.lab.NPP = F}
  if(missing(Spore.influx)){Spore.influx = FALSE}
  if(missing(nlake)){nlake = "Site1"}
  if(missing(NPP.conc.round)){NPP.conc.round =2}
  if(missing(Seuil.NPP)){Seuil.NPP = 300}
  if(missing(Time.window)){Time.window = 500}
  if(missing(Sort.taxon)){Sort.taxon = "Auto"}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  
  
  #### Prepa data ####
  row.names(MP_lake) <- gsub("cf\\.", "cf", row.names(MP_lake))
  row.names(MP_lake) <- gsub("\\.", " ", row.names(MP_lake))
  
  if(is.null(Limites) == F){
    MA_lake <- data.frame(t(MA_lake), check.names = F)
    MA_lake <- MA_lake[MA_lake$Age >= Limites[1] & MA_lake$Age <= Limites[2],]
    MA_lake <- data.frame(t(MA_lake), check.names = F)
    To.keep <- intersect(names(MA_lake), names(MP_lake))
    MP_lake <- MP_lake[To.keep]
  }
  Prof <- MA_lake[1,]
  Prof <- Prof[colnames(MP_lake)]
  Age <- MA_lake["Age",]
  Age <- Age[colnames(MP_lake)]
  
  #### Fusion taxons rares #### 
  MP_max <- MP_lake[0]
  MP_max[, "max"] <- apply(MP_lake[,], 1, max)
  MP_dom <- MP_lake[rowSums(MP_max)>Seuil.NPP,]
  MP_rare <- MP_lake[rowSums(MP_max)<=Seuil.NPP,]
  MP_dom["Other Spores" ,] <- colSums(MP_rare)
  
  #### Couleurs #### 
  InfoPol  <- read.csv(file = Path.index, sep=",",dec=".", header = T, stringsAsFactors = F, row.names = 1)
  if(any(unique(InfoPol$Cortege) %in% c("Foret", "Tourbiere", "PaturageFort", "PaturageMoyen", "Erosion", "Parasitique")) == T){Old.cortege = T}
  else{Old.cortege = F}
  
  if(Col.cortege == F){
    InfoPol  <- cbind(InfoPol, Couleur=rep(InfoPol$AP.NAP))
    InfoPol$Couleur<- gsub('NAP', '#f27041ff', InfoPol$Couleur)
    InfoPol$Couleur<- gsub('AP', '#6db38fff', InfoPol$Couleur)
    InfoPol$Couleur<- gsub('NaN', '#cfce90ff', InfoPol$Couleur)
    InfoPol$Couleur<- gsub('Spore', '#aa373aff', InfoPol$Couleur)
    InfoPol$Couleur<- gsub('Algue', '#62cfdaff', InfoPol$Couleur)
    InfoPol$Couleur<- gsub('Sum', '#aa373aff', InfoPol$Couleur)
  }
  else{
    InfoPol  <- cbind(InfoPol, Couleur=rep(InfoPol$Cortege))
    if(Old.cortege == T){
      InfoPol$Couleur<- gsub('Foret', '#779d3e', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Tourbiere', '#6bbab4', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('PaturageFort', '#ee4a1d', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('PaturageMoyen', '#ee4a1d', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Erosion', '#eeb448ff', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Parasitique', 'grey', InfoPol$Couleur)
    }
    else{
      InfoPol$Couleur[grepl('Saprophytic spores', InfoPol$Couleur)] <- '#779d3e'
      InfoPol$Couleur[grepl('Hydrophytic', InfoPol$Couleur)] <- '#6bbab4'
      InfoPol$Couleur<- gsub('Wetland markers', '#6bbab4', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Dung fungal spores', '#ee4a1d', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Erosion markers', '#ee883f', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Aridity markers', '#eeb448', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Other NPPs', 'grey', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Parasitic spores', 'grey', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Plant pathogens', 'grey', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Ubiquitous spores', 'grey10', InfoPol$Couleur)
    }
    
    InfoPol$Couleur<- gsub('Spore Diversity', '#aa373aff', InfoPol$Couleur)
    InfoPol$Couleur<- gsub('NaN', '#aa373aff', InfoPol$Couleur)
    InfoPol$Couleur[InfoPol$Couleur == ""] <- "grey"
    InfoPol$Couleur[is.na(InfoPol$Couleur)] <- "grey"
  }
  
  if(length(setdiff(row.names(MP_lake),InfoPol$Nom))>0){
    print("**** Attention: these taxa are not in the Index. ****")
    print(setdiff(row.names(MP_lake),InfoPol$Nom))
  }
  
  # print(InfoPol)
  
  #### Calcul diversite ####
  if(Spore.diversity == TRUE){
    if(is.null(Pol.Sum) == T | Diversity.pollen.ratio == F){MatDiv <- Diversite_pollen(MP_lake)}
    if(is.null(Pol.Sum) == F & Diversity.pollen.ratio == T){MatDiv <- Diversite_pollen(MP_lake)
    MatDiv <- 10000*MatDiv/Pol.Sum}}
  else{MatDiv <- NULL}
  
  
  
  #### Calcul Richesse en Type ####
  if(Richesse.type == T ){
    if(Old.cortege == T){
      MForet <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Foret"]),])
      MTourbiere <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Tourbiere"]),])
      MPaturageFort <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "PaturageFort"]),])
      MPaturageMoyen <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "PaturageMoyen"]),])
      MPaturage <- MPaturageMoyen + MPaturageFort
      MErosion <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Erosion"]),])
      MParasite <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Parasitique"]),])
      MatDiv2 <- rbind(Saprophite = MForet, Peat = MTourbiere, Erosion = MErosion, Parasite = MParasite, Grazing = MPaturage)
    }
    else{
      MForet <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Saprophytic spores"]),])
      MTourbiere <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Wetland markers"]),])
      MPaturage <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Dung fungal spores"]),])
      MErosion <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Erosion markers"]),])
      MParasite <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Plant pathogens"]),])
      Maridity <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Aridity markers"]),])
      Mother <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Other NPPs"]),])
      MatDiv2 <- rbind(Saprophilous_Spore = MForet, Peat_Spore = MTourbiere, Erosion_Spore = MErosion, Parasitic_Spore = MParasite, Grazzing_Spore = MPaturage, Aridity_Spore = Maridity, "Other NPPs" = Mother)
    }
    
    if(is.null(Pol.Sum) == F & Diversity.pollen.ratio == T){MatDiv2 <- 10000*MatDiv2/Pol.Sum}
    if(Spore.diversity == T){MatDiv <- rbind(MatDiv, MatDiv2)}
    if(Spore.diversity == F){MatDiv <- MatDiv2}
  }
  
  #### Calcul Cortege ####
  if(Cortege == T){
    if(Old.cortege == T){
      Foret <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Foret"]),])
      Tourbiere <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Tourbiere"]),])
      PaturageFort <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "PaturageFort"]),])
      PaturageMoyen <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "PaturageMoyen"]),])
      Erosion <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Erosion"]),])
      Parasite <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Parasitique"]),])
      PaturageMoyen <-PaturageMoyen + PaturageFort
      Mcortege <- rbind(Saprophilous_Spore = Foret, Peat_Spore = Tourbiere, Erosion_Spore = Erosion, Parasitic_Spore = Parasite, Grazzing_Spore = PaturageMoyen)
      ColCor <- data.frame(Couleur = c("#779d3eff","#6bbab4ff","#eeb448ff","grey", "#ee4a1dff"))
    }  
    else{
      Foret <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Saprophytic spores"]),])
      Tourbiere <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Wetland markers"]),])
      Paturage <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Dung fungal spores"]),])
      Erosion <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Erosion markers"]),])
      Other <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Other NPPs"]),])
      Parasite <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Plant pathogens"]),])
      Aridity <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Aridity markers"]),])
      Mcortege <- rbind("Other NPPs" = Other, Peat_Spore = Tourbiere, Saprophilous_Spore = Foret, Parasitic_Spore = Parasite, Erosion_Spore = Erosion, Grazzing_Spore = Paturage, Aridity_Spore = Aridity)
      ColCor <- data.frame(Couleur = c("#779d3e","#6bbab4","purple","#ee883f","#ee4a1d", "#eeb448", "grey"))
    }  
    row.names(ColCor)<-row.names(Mcortege)
    if(is.null(Keep.cortege) == F){
      Mcortege <- Mcortege[row.names(Mcortege) %in% Keep.cortege,]
      ColCor <- subset(ColCor, row.names(ColCor) %in% Keep.cortege, Couleur)
    }
    
    Keep.taxa.lab <- c(row.names(MP_dom), row.names(Mcortege))
    spec <- cbind(data.frame(t(MP_dom), check.names = F), data.frame(t(Mcortege), check.names = F))
    names(spec) <- Keep.taxa.lab
  }
  if(Cortege == F){
    Keep.taxa.lab <- row.names(MP_dom)
    spec <- cbind(data.frame(t(MP_dom), check.names = F))
    names(spec) <- Keep.taxa.lab
    Mcortege = MP_dom[0]
    ColCor = NULL}
  
  spec <- spec[order(as.double(Age)),]
  Age <- Age[order(as.double(Age))]
  
  #### Calcul des pourcentage ####
  if(Conc.Plot == T){
    if(is.null(Pol.Sum) == T){spec <- (spec*100)/rowSums(MP_lake)}
    if(is.null(Pol.Sum) == F){spec <- (spec*100)/(rowSums(MP_lake)+Pol.Sum)}
  }
  
  #### Sort.Taxon ####
  Sort.col <- function(Type){
    if(missing(Type)){Type = "Normal"}
    if(Type == "Auto.alpha"){InfoPol <- InfoPol[order(InfoPol$Nom),]}
    if(Type == "Auto.mean"){
      spec <- spec[,order(colMeans(spec), decreasing = T)]
      InfoPol <- InfoPol[match(names(spec), InfoPol$Nom),]
    }
    InfoPol <- InfoPol[order(InfoPol$Couleur),]
    inter <- intersect(InfoPol$Nom, names(spec))
    return(inter)}
  
  if(Sort.taxon == 'None'){print("Not sorted")}
  if(Sort.taxon == 'Auto'){
    spec <- spec[,match(Sort.col(),names(spec))]}
  if(Sort.taxon == 'Auto.alpha'){
    spec <- spec[,match(Sort.col(Type = "Auto.alpha"),names(spec))]}    
  if(Sort.taxon == 'Auto.mean'){
    spec <- spec[,match(Sort.col(Type = "Auto.mean"), names(spec))]}
  if(Sort.taxon == 'Manual'){
    if(missing(Manual.sort)){print("**ERROR** The sorting vector is not working. We applied an auto-sorting.")
      spec <- spec[,match(Sort.col(),names(spec))]}
    else{spec <- spec[,match(Manual.sort,names(spec))]}}
  
  #### CONISS Clustering ####
  if(CONISS == TRUE){
    diss <- dist(sqrt(t(MP_dom)/100)^2)
    Clust_DP <- chclust(diss, method = "coniss")
  }
  else{Clust_DP = NULL}
  
  
  
  #### Calcul influx spore ####
  if(Conc.Plot == F & Spore.influx == T){
    Minflux <- rowSums(data.frame(t(MP_lake)))
    Minflux <- Minflux[names(Age)]
    spec["Spore.Influx"] = Minflux
  }
  
  #### Couleur  #### 
  if(Col.cortege == F){
    C <- subset(InfoPol, select=c(Couleur))
    row.names(C)<- InfoPol$Nom
    
    MP_color <- cbind(MP_dom[0], Couleur = C[rownames(MP_dom),])
    MP_color <- rbind(MP_color, ColCor)
    couleur = as.character(MP_color[[1]])
    if(Conc.Plot == F & Spore.influx == T){couleur[length(couleur)+1] = "#aa373aff"}}
  else{couleur <- InfoPol[match(names(spec), InfoPol$Nom), "Couleur"]}
  
  #### Add.plot concentration, NPP ####
  if(is.null(MatDiv) == T){P1L = 0
  P1R = 1}
  
  if(is.null(MatDiv) == F){P1L = 0
  if(length(MatDiv[,1]) == 1){P1R = 0.86}
  else{P1R = (1 - 0.07*length(MatDiv[,1]))}
  P2L = P1R
  P2R = 1}
  
  #### Graph settings (Log trans, auto.lab, Xticks ####
  if(Auto.lab.NPP == T & is.null(Path.index) == F){
    NPP.lab <- InfoPol$Label[match(gsub("\\.", " ", names(spec)), InfoPol$Nom)]
  }
  else{NPP.lab <- NULL}
  
  X.ticks <- round(max(spec), digits = -NPP.conc.round)/4
  if(Log.trans == T){
    spec <- log10(spec*100)
    X.ticks <- round(max(spec), digits = 0)/4
  }
  
  if(is.null(Time.window) == F){Time.window <- seq(-50,round(max(Age),digits=-2),Time.window)}
  
  if(is.null(Order.fungal) == F){print(names(spec))
    spec <- spec[Order.fungal]
    NPP.lab <- NPP.lab[Order.fungal]
    couleur <- couleur[Order.fungal]
  }
  
  #### Plot ####
  if(Conc.Plot == F){
    p <- strat.plot(spec, yvar = as.double(Age), 
                    y.rev = TRUE, 
                    scale.percent = T,     # True pour des pourcentages
                    srt.xlabel = 45,           # Rotation de 45 des noms de taxon
                    xSpace = 0.006,
                    ylabel = "Time (yr cal BP)",
                    y.tks = Time.window,
                    x.pc.inc = X.ticks,
                    plot.poly = FALSE, 
                    plot.bar = TRUE,
                    plot.line = FALSE, 
                    lwd.bar = 7,
                    col.bar = couleur,#      "#8c92a1ff",         # Fait apparaitre les traits continus
                    exag=FALSE,                # Fait apparaitre la zone x10
                    col.exag="auto",           # Zone x10 couleur auto
                    clust = Clust_DP,
                    clust.width=0.075,
                    xRight = P1R,
                    xLeft = P1L
    )}
  if(Conc.Plot == T){
    p <- strat.plot(spec, yvar = as.double(Age), 
                    y.rev = TRUE, 
                    xSpace = 0.002,
                    scale.percent = TRUE,     # True pour des pourcentages
                    srt.xlabel = 45,           # Rotation de 45 des noms de taxon
                    ylabel = "Time (yr cal BP)",
                    y.tks = Time.window,
                    x.pc.inc = round(max(spec, na.rm = F), digits = -3)/4,
                    plot.poly = TRUE, 
                    plot.bar = TRUE,
                    plot.line = FALSE, 
                    lwd.bar = 1,
                    col.poly = couleur,#      "#a2dbe1ff",         # bleu clair
                    col.bar = couleur,#      "#8c92a1ff",         # Fait apparaitre les traits continus
                    col.poly.line = "#8c92a1ff",         # gris fonce
                    exag=T,                # Fait apparaitre la zone x10
                    col.exag="auto",           # Zone x10 couleur auto
                    clust = Clust_DP,
                    clust.width=0.075, 
                    xRight = P1R,
                    xLeft = P1L
    )}
  if(CONISS == T & Nzone > 0){addClustZone(p, Clust_DP, Nzone, 
                                           lwd=1.5, lty=2, col="grey25") 
  }
  if(is.null(MatDiv) == F){
    MatDiv <- data.frame(t(MatDiv))
    MatDiv <- MatDiv[names(Age),]
    names(MatDiv)[1] <- "Spore Diversity"
    
    if(Richesse.type == F){couleur.div = "#aa373aff"}
    if(Richesse.type == T){couleur.div = c("#aa373aff", as.character(ColCor[[1]]))}
    
    strat.plot(MatDiv, yvar = as.double(Age),
               y.rev = TRUE, 
               scale.percent = T,     # True pour des pourcentages
               srt.xlabel = 45,           # Rotation de 45 des noms de taxon
               xSpace = 0.007,
               x.pc.inc = round(max(MatDiv), digits = -1)/2,
               y.axis        = F,
               plot.poly     = F, 
               plot.bar      = F,
               plot.line     = T, 
               lwd.bar       = 1,
               lwd.line = 2,
               col.poly      = couleur.div,          # gris clair
               col.bar       = "#8c92a1ff",         # Fait apparaitre les traits continus
               col.poly.line = couleur.div,     # gris fonce
               col.line = couleur.div,
               exag = F,                # Fait apparaitre la zone x10
               col.exag="auto",           # Zone x10 couleur auto
               xLeft = P2L,
               xRight = P2R,
               add = T
    )}
  if(is.null(Save.plot) == F){dev.off()}
  
  #### Return data ####
  return(list(spec, Age, MatDiv, InfoPol))
}

DP_pol <- function(MP, MA_lake = NULL, Plot.x = NULL, Malg_conc = NULL, Mpol_conc = NULL, Malg_frac, NPP_lake, Keep.pollen = NULL,
                   Path.index = NULL, Save.plot, W, H, Order.algues = NULL, Keep.cortege = NULL, Order.fungal = NULL, Keep.algue = NULL,
                   nlake, Seuil.Pour, Seuil.Pour.alg, Seuil.Conc.alg, Seuil.NPP, Auto.Pollen.lab, Order.pollen = NULL, Keep.NPP = NULL,
                   Sort.taxon, Pollen.lab, Cortege, Richesse.type, Only.AP, Zone.hiatus, Time.window, Show.stats = T,
                   Save.missing.tax = NULL,
                   Manual.sort = NULL, AP.NAP, CONISS, Nzone, Pollen.diversity, Spore.diversity, Limites, Al.Pol.influx,
                   Pol.influx.total, Spore.influx, Ratio, Save.path, Zone.clim, Name.zone, Temp.zone){
  #### Default values ####
  if(missing(Limites)){Limites = NULL}
  if(missing(Zone.clim)){Zone.clim = NULL}
  if(missing(Name.zone)){Name.zone = NULL}
  if(missing(Temp.zone)){Temp.zone = rep("C", length(Zone.clim))}
  if(missing(Cortege)){Cortege = FALSE}
  if(missing(Richesse.type)){Richesse.type = FALSE}
  if(missing(Malg_conc) & missing(Malg_frac)){Plot.algue = F}
  else{Plot.algue = T}
  if(missing(Malg_frac)){Malg_frac = NULL}
  
  if(missing(NPP_lake)){Plot.NPP = F}else{Plot.NPP = TRUE}
  if(missing(AP.NAP)){AP.NAP = FALSE}
  if(missing(Only.AP)){Only.AP = FALSE}
  if(missing(CONISS)){CONISS = FALSE}
  if(missing(Nzone)){Nzone = FALSE}
  Diversite.plot = F
  if(missing(Pollen.diversity)){Pollen.diversity = FALSE}
  if(missing(Spore.diversity)){Spore.diversity = FALSE}
  if(Pollen.diversity == T | Spore.diversity == T){Diversite.plot = T}
  if(missing(Pol.influx.total)){Pol.influx.total = FALSE}
  if(missing(Al.Pol.influx)){Al.Pol.influx = FALSE}
  if(missing(Spore.influx)){Spore.influx = FALSE}
  if(missing(nlake)){nlake = NULL}
  if(missing(Seuil.Pour)){Seuil.Pour = 0.01}
  if(missing(Seuil.Pour.alg)){Seuil.Pour.alg = 0.01}
  if(missing(Seuil.Conc.alg)){Seuil.Conc.alg = 2000}
  if(missing(Seuil.NPP)){Seuil.NPP = 300}
  if(missing(Time.window)){Time.window = NULL}
  if(missing(Sort.taxon)){Sort.taxon = "Auto"}
  NBDisplay <- c(Plot.algue, Plot.NPP, Pollen.diversity, Pol.influx.total)
  if(missing(Ratio)){Ratio = NULL}
  if(missing(Zone.hiatus)){Zone.hiatus = NULL}
  if(missing(Pollen.lab)){Pollen.lab = NULL}
  if(missing(Auto.Pollen.lab)){Auto.Pollen.lab = F}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  
  #### Graphical settings ####
  if(is.null(Zone.clim) ==F){
    Rect.zone = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)],
                           xmax = Zone.clim[seq(2,length(Zone.clim), by=2)],
                           Temp.col = Temp.zone)}
  else{Rect.zone = data.frame(xmin = 0, xmax = 0, Temp.col = "black")}
  
  if(is.null(Name.zone) == F){Rect.zone = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)],
                                                     xmax = Zone.clim[seq(2,length(Zone.clim), by=2)],
                                                     Temp.col = Temp.zone,
                                                     Title.zone = Name.zone)}
  else{Rect.zone = data.frame(xmin = 0, xmax = 0, Temp.col = "", Title.zone = "")}
  
  Rect.zone$Temp.col <- gsub('C', '#E5EDFB', Rect.zone$Temp.col)
  A <- gsub('#E5EDFB', '#336CAD', Rect.zone$Temp.col)
  Rect.zone$Temp.col <- gsub('W', '#F8E2E3', Rect.zone$Temp.col)
  A <- gsub('W', '#D55C66', A)
  Rect.zone$Temp.col <- gsub('G', '#BEBEBE', Rect.zone$Temp.col)
  A <- gsub('G', 'grey10', A)
  Rect.zone$Name.col <- A
  
  #### MA (or Plot x) missing ####
  All.ages <- c("Age", "agebp", "MAge", "AgeBP2023")
  All.depths <- c("Depth", "Top", "Bottom", "Prof", "Profondeur")
  
  if(is.null(MA_lake) == T){
    print("No age-depth model provided, nor ordination. Samples ploted by alphabetical order.")
    MA_lake <- setNames(data.frame(t(data.frame(Age = seq(1:ncol(MP))))), names(MP))
    Plot.x <- "Ordin"
  }
  else{
    if(is.null(Plot.x) == T){
      Plot.x <- c(All.ages, All.depths)
      Plot.x <- row.names(MA_lake)[match(Plot.x, row.names(MA_lake))[1]]
      
      
      if(is.na(Plot.x) == T){
        print("**** No automatic plot x variable found in MA_lake. Please provide one (row.names(MA_lake)). ****")
        MA_lake <- setNames(data.frame(t(data.frame(Age = seq(1:ncol(MP))))), names(MP))
      }
      else{print(paste("****", Plot.x, "selected by default for plot x variable. ****"))}
    }
  }
  
  #### Prepa data ####
  Prep.data <-function(MP_lake, Max_seuil, Type, Calc.stat = F, Only.keep = NULL, Save.missing.tax = NULL){
    #### Default settings ####
    if(missing(Type)){Type = "Unknown"}
    
    #### Limites ####
    if(is.null(Limites) == F){
      MA_lake <- data.frame(t(MA_lake))
      MA_lake <- MA_lake[MA_lake$Age >= Limites[1] & MA_lake$Age <= Limites[2],]
      MA_lake <- data.frame(t(MA_lake))
      To.keep <- intersect(names(MA_lake), names(MP_lake))
      MP_lake <- MP_lake[To.keep]
    }
    
    #### Init Matrice + Age selection ####
    Age <- MA_lake[Plot.x,]
    Age <- Age[colnames(MP_lake)]
    Age <- Age[order(as.double(Age))]
    MP_lake <- MP_lake[colnames(Age)]
    
    #### Couleurs #### 
    if(is.null(Path.index) == F){
      InfoPol  <- read.csv(file = Path.index, sep=",",dec=".", header = T, stringsAsFactors = F, row.names = 1, check.names = FALSE)
      InfoPol  <- cbind(InfoPol, Couleur=rep(InfoPol$AP.NAP))
      InfoPol$Couleur<- gsub('NAP', '#f27041ff', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('AP', '#6db38fff', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('NaN', '#cfce90ff', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Unknown', '#cfce90ff', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Spore', '#aa373aff', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Algue', '#4666E9', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Algal', '#4666E9', InfoPol$Couleur)
      InfoPol$Couleur<- gsub('Sum', '#aa373aff', InfoPol$Couleur)
      InfoPol$Couleur[grep('Shrub', InfoPol$Cortege)] <- '#e0b300ff'
      
      if(length(setdiff(row.names(MP_lake),InfoPol$Nom))>0){
        print("**** Attention: these taxa are not in the Index. ****")
        print(setdiff(row.names(MP_lake),InfoPol$Nom))
        if(is.null(Save.missing.tax) == F){
          write.table(setdiff(row.names(MP_lake),InfoPol$Nom), Save.missing.tax)}
      }
    }
    else{
      InfoPol <- data.frame(Index = seq(1,nrow(MP_lake)),
                            Nom = row.names(MP_lake),
                            Label = row.names(MP_lake),
                            Type = "Pollen", AP.NAP = NaN, Cortege = NaN, Couleur = "#6db38fff"
      )
    }
    
    #### Calcul diversite ####
    if(Pollen.diversity == TRUE | Spore.diversity == TRUE){MatDiv <- Diversite_pollen(MP_lake)}
    else{MatDiv <- NULL}
    
    #### Special Fungi ####
    if(Type == "NPP"){
      if(any(unique(InfoPol$Cortege) %in% c("Foret", "Tourbiere", "PaturageFort", "PaturageMoyen", "Erosion", "Parasitique")) == T){Old.cortege = T}
      else{Old.cortege = F}
      
      #### Calcul Richesse en Type ####
      if(Richesse.type == T ){
        if(Old.cortege == T){
          MForet <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Foret"]),])
          MTourbiere <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Tourbiere"]),])
          MPaturageFort <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "PaturageFort"]),])
          MPaturageMoyen <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "PaturageMoyen"]),])
          MPaturage <- MPaturageMoyen + MPaturageFort
          MErosion <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Erosion"]),])
          MParasite <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Parasitique"]),])
          MatDiv2 <- rbind(Saprophite = MForet, Peat = MTourbiere, Erosion = MErosion, Parasite = MParasite, Grazing = MPaturage)
        }
        else{
          MForet <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Saprophytic spores"]),])
          MTourbiere <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Wetland markers"]),])
          MPaturage <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Dung fungal spores"]),])
          MErosion <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Erosion markers"]),])
          MParasite <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Parasitic spores"]),])
          Maridity <- Diversite_pollen(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Aridity markers"]),])
          MatDiv2 <- rbind(Saprophite = MForet, Peat = MTourbiere, Erosion = MErosion, Parasite = MParasite, Grazing = MPaturage, Aridity = Maridity)
        }
        
        if(Spore.diversity == T){MatDiv <- rbind(MatDiv, MatDiv2)}
        if(Spore.diversity == F){MatDiv <- MatDiv2}
      }
      
      #### Calcul Cortege ####
      if(Cortege == T){
        if(Old.cortege == T){
          Foret <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Foret"]),])
          Tourbiere <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Tourbiere"]),])
          PaturageFort <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "PaturageFort"]),])
          PaturageMoyen <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "PaturageMoyen"]),])
          Erosion <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Erosion"]),])
          Parasite <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Parasitique"]),])
          PaturageMoyen <-PaturageMoyen + PaturageFort
          Mcortege <- rbind(Saprophilous_Spore = Foret, Peat_Spore = Tourbiere, Erosion_Spore = Erosion, Parasitic_Spore = Parasite, Grazzing_Spore = PaturageMoyen)
          ColCor <- data.frame(Couleur = c("#779d3eff","#6bbab4ff","#eeb448ff","grey", "#ee4a1dff"))
        }  
        else{
          Foret <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Saprophytic spores"]),])
          Tourbiere <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Wetland markers"]),])
          Paturage <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Dung fungal spores"]),])
          Erosion <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Erosion markers"]),])
          Parasite <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Parasitic spores"]),])
          Aridity <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$Cortege == "Aridity markers"]),])
          Mcortege <- rbind(Saprophilous_Spore = Foret, Peat_Spore = Tourbiere, Erosion_Spore = Erosion, Parasitic_Spore = Parasite, Grazzing_Spore = Paturage, Aridity_Spore = Aridity)
          ColCor <- data.frame(Couleur = c("#779d3eff","#6bbab4ff","purple","#ee883fff","#ee4a1dff", "#eeb448ff"))
        }  
        
        row.names(ColCor)<-row.names(Mcortege)
        
        if(is.null(Keep.cortege) == F){
          Mcortege <- Mcortege[row.names(Mcortege) %in% Keep.cortege,]
          ColCor <- subset(ColCor, row.names(ColCor) %in% Keep.cortege, Couleur)
        }
        
      }
      
      if(Cortege == F){
        Mcortege = MP_dom[0]
        ColCor = NULL}
    }
    
    #### Fusion taxons rares ####
    MP_max <- MP_lake[0]
    MP_max[, "max"] <- apply(MP_lake[,], 1, max)
    MP_dom <- MP_lake[rowSums(MP_max)>Max_seuil,]
    MP_rare <- MP_lake[rowSums(MP_max)<=Max_seuil,]
    
    if(is.null(Only.keep) == F){
      Match.OK <- Only.keep[which(Only.keep %in% row.names(MP_lake))]
      Wrong <- Only.keep[which(!Only.keep %in% row.names(MP_lake))]
      if(length(Wrong) > 0){
        print("****"); print("The following taxa are wrong in the Keep.only vector:")
        print(Wrong)
        print("****")}
      
      MP_dom <- MP_lake[row.names(MP_lake) %in% Match.OK,]
      MP_rare <- MP_lake[!row.names(MP_lake) %in% Match.OK,]
    }
    
    AP.rare <- MP_rare[intersect(row.names(MP_rare),InfoPol$Nom[InfoPol$AP.NAP == "AP"]),]
    NAP.rare <- MP_rare[intersect(row.names(MP_rare),InfoPol$Nom[InfoPol$AP.NAP == "NAP"]),]
    
    if(length(row.names(AP.rare)) >= 1){
      MP_dom["Other Trees" ,] <- colSums(AP.rare)
      MP_dom["Other Herbs" ,] <- colSums(NAP.rare)
    }
    else{MP_dom["Other" ,] <- colSums(MP_rare)}
    
    #### Sort.Taxon ####
    Sort.col <- function(Type){
      if(missing(Type)){Type = "Normal"}
      if(Type == "Auto.alpha"){InfoPol <- InfoPol[order(InfoPol$Nom),]}
      if(Type == "Auto.mean"){
        MP_dom <- MP_dom[order(rowMeans(MP_dom), decreasing = T),]
        InfoPol <- InfoPol[match(row.names(MP_dom), InfoPol$Nom),]
      }
      InfoPol <- InfoPol[order(InfoPol$Couleur),]
      inter <- intersect(InfoPol$Nom, row.names(MP_dom))
      return(inter)}
    
    if(Sort.taxon == 'None'){print("Not sorted")}
    if(Sort.taxon == 'Auto'){
      MP_dom <- MP_dom[match(Sort.col(),row.names(MP_dom)),]}
    if(Sort.taxon == 'Auto.alpha'){
      MP_dom <- MP_dom[match(Sort.col(Type = "Auto.alpha"),row.names(MP_dom)),]}    
    if(Sort.taxon == 'Auto.mean'){
      MP_dom <- MP_dom[match(Sort.col(Type = "Auto.mean"),row.names(MP_dom)),]}
    if(Sort.taxon == 'Manual'){
      if(is.null(Manual.sort) == T){print("**ERROR** The sorting vector is not working. We applied an auto-sorting.")
        MP_dom <- MP_dom[match(Sort.col(),row.names(MP_dom)),]}
      else{MP_dom <- MP_dom[match(Manual.sort,row.names(MP_dom)),]}}
    
    #### Calcul AP et NAP  #### 
    if(AP.NAP == TRUE){
      AP <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$AP.NAP == "AP"]),])
      NAP <- colSums(MP_lake[intersect(row.names(MP_lake),InfoPol$Nom[InfoPol$AP.NAP == "NAP"]),])
      AP_nul <- sum(AP)
      if(AP_nul != 0){
        if(Only.AP == F){MP_dom <- rbind(AP = AP, NAP = NAP, MP_dom)}
        else{MP_dom <- rbind(AP = AP, MP_dom)}
      }
    }  
    
    #### Couleur  #### 
    C <- subset(InfoPol, select=c(Couleur))
    row.names(C)<- InfoPol$Nom
    MP_color <- cbind(MP_dom[0], Couleur = C[rownames(MP_dom),])
    if(Type == "NPP"){
      MP_color <- rbind(MP_color, ColCor)}
    couleur = as.character(MP_color[[1]])
    
    #### CONISS Clustering #### 
    if(CONISS == TRUE){
      diss <- dist(sqrt(t(MP_dom)/100)^2)
      Clust_DP <- chclust(diss, method = "coniss")
      
      if(Calc.stat == T){
        cc <- cutree(Clust_DP, k = Nzone)
        Tab.stat <- data.frame(t(rbind(CONISS = cc, MP_lake)))
        
        Age.CONISS <- data.frame(cbind(Age = t(Age), CONISS = Tab.stat$CONISS))
        Age.CONISS.min <- aggregate(Age.CONISS,  by = Age.CONISS["CONISS"], FUN = min)
        Age.CONISS.max <- aggregate(Age.CONISS,  by = Age.CONISS["CONISS"], FUN = max)
        Age.CONISS <- cbind(Age.CONISS.min[c(1:2)], Age.CONISS.max[c(2)])
        Age.CONISS$CONISS <- paste("U", Age.CONISS$CONISS, sep = "_") 
        names(Age.CONISS) <- c("CONISS", "Age_min", "Age_max")
        
        
        Tab.stat <- aggregate(Tab.stat, list(Tab.stat[["CONISS"]]), FUN = mean, na.action = na.omit)
        Tab.stat <- Tab.stat[-c(1)]
        row.names(Tab.stat) <- paste("U", Tab.stat$CONISS, sep = "_")
        Tab.stat <- data.frame(t(round(Tab.stat[-c(1)], digits = 2)*100))
        
        Resume <- c()
        for(i in 1:ncol(Tab.stat)){
          NN = row.names(Tab.stat[Tab.stat[i] > 2,])
          VV = Tab.stat[Tab.stat[i] > 2,i]
          NN = NN[order(VV, decreasing = T)]
          VV = VV[order(VV, decreasing = T)]
          JJ = paste(NN, " (", VV, " %)", sep = '', collapse = ', ')
          Resume[i] <- JJ
        }
        Tab.stat <- data.frame(t(Tab.stat))
        Tab.stat <- cbind(Tab.stat, Resume, Age.CONISS[c(2,3)])
        if(Show.stats == T){print(Resume)}
      }
      else{Tab.stat <- NULL}
    }
    else{
      Tab.stat <- NULL
      Clust_DP = NULL}
    
    #### Calcul des Ratios ####
    
    if(is.null(Ratio) == F & Type == "Pollen"){
      print("Ratio calcul is beginning.")
      MP_cal <- data.frame(t(MP_lake), check.names = F)
      MR <- MP_cal[0]
      for(i in 1:length(TabRa)){
        TaxonDeno <- c(as.character(unique(TabRa[[i]]$Deno)))
        ColDeno <- subset(MP_cal, select = TaxonDeno)
        TaxonNomi <- c(as.character(unique(TabRa[[i]]$Nomi)))
        ColNomi <- subset(MP_cal, select = TaxonNomi)
        MR[i] <- rowSums(ColDeno)/rowSums(ColNomi)
        InitLDeno <- paste(substring(TaxonDeno, 1, 2), collapse = "+")
        InitLNomi <- paste(substring(TaxonNomi, 1, 2), collapse = "+")
        colnames(MR)[i] <- paste(InitLDeno,InitLNomi,sep = "/")
      }
      MR <- MR[order(as.double(Age)),]
    }
    else{MR = NULL}
    
    #### Save to csv ####
    if(is.null(Save.path) == F){
      Path.to.create <- gsub("(.*/).*\\.csv.*","\\1",Save.path)
      dir.create(file.path(Path.to.create), showWarnings = FALSE)
      
      if(is.null(MR) == F){
        MR.save <- cbind(Age = t(Age[order(as.double(Age))]), MR)
        Save.path.ratio <- gsub("\\.csv", "_Ratio.csv", Save.path)
        Save.path.ratio.RDS <- gsub("\\.csv", "_Ratio.Rds", Save.path)
        saveRDS(MR.save, Save.path.ratio.RDS)
        write.table(MR.save, file=Save.path.ratio, row.names=T, col.names=NA, sep=",", dec = ".")
      }
      
      if(Diversite.plot == T){
        if(Pollen.diversity == T){
          Save.path.div <- gsub("\\.csv", "_PolDiversity.csv", Save.path)
          write.table(MatDiv, file=Save.path.div, row.names=T, col.names=NA, sep=",", dec = ".")}
        if(Spore.diversity == T){
          Save.path.div <- gsub("\\.csv", "_SporeDiversity.csv", Save.path)
          write.table(MatDiv, file=Save.path.div, row.names=T, col.names=NA, sep=",", dec = ".")}
      }
      
      if(CONISS == T & Calc.stat == T){
        if(Type == "Pollen"){
          Save.path.div <- gsub("\\.csv", "_statPol_CONISS.csv", Save.path)
          write.table(Tab.stat, file=Save.path.div, row.names=T, col.names=NA, sep=",", dec = ".")}
        if(Type == "Algue"){
          Save.path.div <- gsub("\\.csv", "_statAlg_CONISS.csv", Save.path)
          write.table(Tab.stat, file=Save.path.div, row.names=T, col.names=NA, sep=",", dec = ".")}
        
      }
    }
    
    #### Export data ####
    if(Type == "NPP" & Cortege == T){
      Keep.taxa.lab <- c(row.names(MP_dom), row.names(Mcortege))
      spec <- cbind(data.frame(t(MP_dom)), data.frame(t(Mcortege)))
      names(spec) <- Keep.taxa.lab
    }
    else{
      Keep.taxa.lab <- row.names(MP_dom)
      spec <- data.frame(t(MP_dom))
      names(spec) <- Keep.taxa.lab
    }
    
    spec <- spec*100                            # On passe en pourcentage
    spec <- spec[order(as.double(Age)),]
    
    if(!Plot.x %in% c(All.depths, All.ages)){Age <- MA_lake}
    
    full <- list(spec, Age, couleur, Clust_DP, MatDiv, MR, Tab.stat)  # export function
    return(full)
  }
  
  #### Distance entre graph + Legend en bas ####
  if(Plot.algue == T & Plot.NPP == T & Pol.influx.total == T & Diversite.plot == T & is.null(Ratio) == T){
    RPollen = LInf <- 0.5
    RInf = LAlgue <- 0.52
    RAlgue = LNPP <- 0.7
    RNPP = LDiv <- 0.899
    RDiv = .97
    Legend.haut <- c("Climat", "Pollen", "Algae", "Fungal Spore")
    Ajust.haut <- c(0, LAlgue/2, LAlgue, LNPP)
    Legend.bas <- c("%TP", "%TP+NPP", "10^3.#grains.cm-3", "#Taxa")
    Ajust.bas <- c(0, LAlgue, LNPP, LDiv)}
  
  if(Plot.algue == T & Plot.NPP == T & Pol.influx.total == T & Diversite.plot == T & is.null(Ratio) == F){
    RPollen = LInf <- 0.42
    RInf = LAlgue <- 0.45
    RAlgue = LNPP <- 0.7
    RNPP = LDiv <- 0.85
    RDiv = LRatio <- 0.89
    Legend.haut <- c("Climat", "Pollen", "Algae", "Fungal Spore")
    Ajust.haut <- c(0, LAlgue/2, LAlgue, LNPP)
    Legend.bas <- c("%TP", "%TP+NPP", "10^3.#grains.cm-3", "#Taxa")
    Ajust.bas <- c(0, LAlgue, LNPP, LDiv)}
  if(Plot.algue == T & Plot.NPP == T & Pol.influx.total == F & Diversite.plot == T & is.null(Ratio) == F){
    RPollen = LInf <- 0.43
    RInf = LAlgue <- 0.43
    RAlgue = LNPP <- 0.65
    RNPP = LDiv <- 0.85
    RDiv = LRatio <- 0.89
    Legend.bas <- c("%TP", expression(10^6 ~ "#grain." ~ cm^-3 ~ "& %TP+NPP"), expression("#grain." ~ cm^-3), "#Taxa", "Arbitrary Unit")
    Legend.haut <- c("Climat", "Pollen", "Algae", "Fungal Spore", "Pollen Ratio")
    Ajust.haut <- c(0, LAlgue/2, LAlgue, LNPP, LRatio)
    Ajust.bas <- c(0, LAlgue, LNPP, LDiv, LRatio)}
  if(Plot.algue == T & Plot.NPP == T & Pol.influx.total == F & Diversite.plot == T & is.null(Ratio) == T){
    RPollen = LInf <- 0.47
    RInf = LAlgue <- 0.47
    RAlgue = LNPP <- 0.73
    RNPP = LDiv <- 0.89
    RDiv = LRatio <- 0.99
    Legend.haut <- c("Climat", "Pollen", "Algae", "Fungal Spore")
    Ajust.haut <- c(0, LAlgue/2, LAlgue, LNPP)
    Legend.bas <- c("%TP", "%TP+NPP", "#grains.cm-3", "#Taxa")
    Ajust.bas <- c(0, LAlgue, LNPP, LDiv)}
  if(Plot.algue == T & Plot.NPP == T & Pol.influx.total == F & Diversite.plot == F & is.null(Ratio) == T){
    RPollen = LInf <- 0.55
    RInf = LAlgue <- 0.55
    RAlgue = LNPP <- 0.8
    RNPP = LDiv <- 0.99
    RDiv = LRatio <- 0.99
    Legend.bas <- c("%TP", "%TP+NPP", "#grain.cm-3")
    Legend.haut <- c("Climat", "Pollen", "Algae", "Fungal Spore")
    Ajust.haut <- c(0, LAlgue/2, LAlgue, LNPP)
    Ajust.bas <- c(0, LAlgue, LNPP)} 
  if(Plot.algue == T & Plot.NPP == F & Pol.influx.total == T & Diversite.plot == F & is.null(Ratio) == T){
    RPollen = LInf <- 0.55
    RInf = LAlgue <- 0.6
    RAlgue = LNPP <- 0.8
    RNPP = LDiv <- 0.99
    RDiv = LRatio <- 0.99
    Legend.bas <- c("%TP", "%TP+NPP", "#grain.cm-3")
    Legend.haut <- c("Climat", "Pollen", "Algae", "Fungal Spore", "Pollen Ratio")
    Ajust.haut <- c(0, LAlgue/2, LAlgue, LNPP, LRatio)
    Ajust.bas <- c(0, LAlgue, LNPP)}
  if(Plot.algue == F & Plot.NPP == F & Pol.influx.total == F & Diversite.plot == T & is.null(Ratio) == F){
    RPollen = LInf <- 0.85
    RInf = LAlgue <- 0.85
    RAlgue = LNPP <- 0.85
    RNPP = LDiv <- 0.85
    RDiv = LRatio <- 0.89
    Legend.bas <- c("%TP", "%TP+NPP", "#grain.cm-3")
    Legend.haut <- c("Climat", "Pollen", "Pollen Ratio")
    Ajust.haut <- c(0, LAlgue/2, LAlgue, LNPP, LRatio)
    Ajust.bas <- c(0, LAlgue, LNPP)
  }
  if(Plot.algue == F & Plot.NPP == F & Pol.influx.total == F & Diversite.plot == T & is.null(Ratio) == T){
    RPollen = LInf <- 0.89
    RInf = LAlgue <- 0.8
    RAlgue = LNPP <- 0.8
    RNPP = LDiv <- 0.89
    RDiv = LRatio <- 0.95
    Legend.bas <- c("%TP")
    Legend.haut <- c("Pollen", "CONISS")
    Ajust.haut <- c(0, LRatio)
    Ajust.bas <- c(0)
  }
  if(Plot.algue == F & Plot.NPP == F & Pol.influx.total == T & Diversite.plot == F & is.null(Ratio) == T){
    RPollen = LInf <- 0.89
    RInf = 0.99
    Legend.bas <- c("%TP", "%TP+NPP", "#grain.cm-3")
    Legend.haut <- c("Climat", "Pollen", "Pollen Ratio")
  }
  if(Plot.algue == F & Plot.NPP == F & Pol.influx.total == T & Diversite.plot == T & is.null(Ratio) == T){
    RPollen = LInf <- 0.85
    RInf = LDiv <- 0.89
    RDiv = LRatio <- 0.95
    Legend.bas <- c("%TP", "%TP+NPP", "#grain.cm-3")
    Legend.haut <- c("Climat", "Pollen", "Pollen Ratio")
  }
  if(Plot.algue == F & Plot.NPP == F & Pol.influx.total == F & Diversite.plot == F & is.null(Ratio) == T){RPollen = LInf <- 1}
  
  #### Mono / double plot ####
  Keep.all.x <- Plot.x
  for(i in Keep.all.x){
    Plot.x <- i 
    #### Save plots ####
    if(is.null(Save.plot) == F){
      Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
      dir.create(file.path(Path.to.create), showWarnings = FALSE)
      if(length(Keep.all.x) > 1){
        print(paste("Plot du DP in", Plot.x))
        Save.plot.i <- paste(gsub("\\.pdf", "", Save.plot), "_", Plot.x, ".pdf", sep = "")
        }
      else{Save.plot.i <- Save.plot}
      
      if(is.null(W) == F & is.null(H) == F){
        pdf(file = Save.plot.i, width = W*0.01041666666667, height = H*0.01041666666667)}
      else{pdf(file = Save.plot.i)}}
    
    #### Plot  ####
    yo <- Prep.data(MP, Seuil.Pour, Only.keep = Keep.pollen, "Pollen", Calc.stat = T, Save.missing.tax = Save.missing.tax)
    spec = yo[[1]]
    Age = yo[[2]]
    couleur = yo[[3]]
    Clust_DP = yo[[4]]
    Ratio = yo[[6]]
    if(Diversite.plot == T | Plot.algue == T | Plot.NPP == T | Pol.influx.total == T){Clust_DP_plot = NULL}
    else(Clust_DP_plot <- Clust_DP)
    
    #### Main plot graph settings ####
    if(Auto.Pollen.lab == T & is.null(Pollen.lab) == T & is.null(Path.index) == F){
      InfoPol  <- read.csv(file = Path.index, sep=",",dec=".", header = T, stringsAsFactors = F, row.names = 1, check.names = FALSE)
      Pollen.lab <- InfoPol$Label[match(gsub("\\.", " ", names(spec)), InfoPol$Nom)]
    }
    
    if(Plot.x %in% All.ages){
      if(is.null(Time.window) == T){Time.window <- 250}
      My_ticks <- seq(-50,round(max(Age),digits=-2),Time.window)
      My_y_lab <- "Time (yr cal BP)"
    }
    if(Plot.x %in% All.depths){
      if(is.null(Time.window) == T){Time.window <- 10}
      My_ticks <- seq(round(min(Age),digits=-1),round(max(Age),digits=-1),Time.window)
      My_y_lab <- "Depth (cm)"
    }
    if(!Plot.x %in% c(All.depths, All.ages)){
      My_ticks <- seq(1:length(Age))
      My_y_lab <- "Sample ordination (#)"
    }
    
    #### Age Zone ####
    Display.age.rioja = F
    Lab.rect.x1 = -0.03
    Lab.rect.x2 = -0.04
    par(mar=c(2.3, #bottom
              3, #left 0 parfait pour plus tard
              9.92, #top
              0)) # rigth
    
    p2 <- plot(10,0,
               bty = "n",
               ylab = My_y_lab,
               xlab = "",
               ylim = c(max(Age, na.rm = T),min(Age, na.rm = T)), yaxt = "n",
               xlim = c(0,1), xaxt = "n",
               yaxp = c(round(max(Age, na.rm = T),digits=-2),
                        round(min(Age, na.rm = T),digits=-2), 20)
    )
    
    Zone.age = T
    if(nrow(Rect.zone) > 1 & is.null(Rect.zone$Temp.col) == F){
      if(is.null(nlake) == F){
        title(main = paste("Pollen diagram of", nlake, "Lake"), line = 7, cex.main = 2.7, font.main = 1, adj = 0)
        My_lake_name <- paste("Pollen diagram of", nlake, "Lake")
      }
      else{My_lake_name <- NULL}
      
      axis(2, at = seq(-50,round(max(Age),digits=-2),Time.window), las = 1, cex.axis = 0.8)
      
      for(z in 1:nrow(Rect.zone)){
        rect(xleft = Lab.rect.x1, xright = 1, ybottom = Rect.zone[z,1], ytop = Rect.zone[z,2], col = Rect.zone[z,3], lty = 3)
        rect(xleft = Lab.rect.x1, xright = Lab.rect.x2, ybottom = Rect.zone[z,1], ytop = Rect.zone[z,2], col = Rect.zone[z,5],)
        text((Lab.rect.x1+Lab.rect.x2)/2, (Rect.zone[z,1] + Rect.zone[z,2])/2, font = 2,
             labels = Rect.zone[z,4], col = "white", srt = 90, cex = 1)
      }
      
      mtext(side = 1, line = 1, text = Legend.bas, adj = Ajust.bas, cex = 1.4)
      mtext(side = 3, line = 4.8, text = Legend.haut, adj = Ajust.haut, cex = 1.9)
    }
    else{
      Display.age.rioja = T
      if(is.null(nlake) == F){
        My_lake_name <- paste("Pollen diagram of", nlake, "Lake")
      }
      else{My_lake_name <- NULL}
    }
    
    #### Zone hiatus ####
    if(is.null(Zone.hiatus) == F){
      Rect.zone.hiatus = data.frame(xmin = Zone.hiatus[1], 
                                    xmax = Zone.hiatus[2],
                                    Temp.col = "grey30",
                                    Title.zone = "Hiatus", 
                                    Name.col = "grey30")
      
      rect(xleft = Lab.rect.x1, xright = 1, ybottom = Rect.zone.hiatus[1,1], ytop = Rect.zone.hiatus[1,2], col = Rect.zone.hiatus[1,3], lty = 3)
      rect(xleft = Lab.rect.x1, xright = Lab.rect.x2, ybottom = Rect.zone.hiatus[1,1], ytop = Rect.zone.hiatus[1,2], col = Rect.zone.hiatus[1,5],)
      text((Lab.rect.x1+Lab.rect.x2)/2, (Rect.zone.hiatus[1,1] + Rect.zone.hiatus[1,2])/2, font = 2,
           labels = Rect.zone.hiatus[1,4], col = "white", srt = 90, cex = 1)}
    
    #### Order pollen columns ####
    if(is.null(Keep.pollen) == F & is.null(Order.pollen) == T){
      Match.OK <- Keep.pollen[which(Keep.pollen %in% names(spec))]
      if(AP.NAP == T){Match.OK <- c("AP", "NAP", Match.OK)}
      Add <- setdiff(names(spec), Match.OK)
      Save.order <- match(Add, names(spec))
      if(length(Add)>0){Match.OK <- c(Match.OK, Add)}
      Order.pollen <- match(Match.OK, names(spec))
    }
    
    if(is.null(Order.pollen) == F){
      spec <- spec[Order.pollen]
      couleur <- couleur[Order.pollen]
      if(is.null(Pollen.lab) == F){Pollen.lab <- Pollen.lab[Order.pollen]}
    }
    
    #### Plot POLLEN % ####
    p <- strat.plot(spec, yvar = as.double(Age), 
                    y.rev = TRUE, 
                    scale.percent = TRUE,     # True pour des pourcentages
                    srt.xlabel = 45,           # Rotation de 45 des noms de taxon
                    title = My_lake_name,
                    ylabel = My_y_lab,
                    y.tks = My_ticks,
                    xSpace = 0.003,
                    y.axis = Display.age.rioja,
                    x.names = Pollen.lab,
                    mgp = c(0, 0.5, 0.3),
                    add = T,
                    plot.poly     = TRUE,
                    plot.bar      = TRUE,
                    plot.line     = TRUE,
                    lwd.bar       = 1,
                    lwd.poly      = 1,
                    col.poly      = couleur,         # bleu clair
                    #col.bar       = "#a2dbe1ff",         # Fait apparaitre les traits continus
                    col.poly.line = "#8c92a1ff",         # gris fonce
                    exag=TRUE,                # Fait apparaitre la zone x10
                    # exag=ex,   ici on peut mettre une exageration sur seulement les taxa minoritaires
                    exag.alpha=0.30, exag.mult=2.5,
                    col.exag="auto",           # Zone x10 couleur auto
                    clust = Clust_DP_plot,
                    clust.width=0.015, 
                    #fun1=sm.fun                # application ou non du lissage
                    xLeft = 0.03,
                    xRight = RPollen
    )
    if(CONISS == TRUE & Nzone > 0){
      addClustZone(p, Clust_DP, Nzone, lwd=1.8, lty=2, col="grey25")}
    
    #### Plot Ratio pollen ####
    if(is.null(Ratio)==F){
      Clust_DP_plot = Clust_DP
      P.ratio <- strat.plot(Ratio, yvar = as.double(Age),
                            y.rev = TRUE, 
                            scale.percent = F,     # True pour des pourcentages
                            srt.xlabel = 45,           # Rotation de 45 des noms de taxon
                            #xlim = c(20,30),
                            xSpace = 0.002,
                            y.axis        = F,
                            plot.poly     = T, 
                            plot.bar      = F,
                            plot.line     = T, 
                            lwd.bar       = 1,
                            mgp = c(0, 0.5, 0.3),
                            col.poly      = "#E6E6C6",          # gris clair
                            col.bar       = "#8c92a1ff",         # Fait apparaitre les traits continus
                            col.poly.line = "#8c92a1ff",     # gris fonce
                            exag = T,                # Fait apparaitre la zone x10
                            col.exag="auto",           # Zone x10 couleur auto
                            clust = Clust_DP_plot,
                            exag.alpha=0.30, exag.mult=2.5,
                            clust.width = 0.015, 
                            #fun1=sm.fun                # application ou non du lissage
                            xLeft = LRatio,
                            xRight = 0.99,
                            add = T
      )
      if(CONISS == TRUE & Nzone > 0){addClustZone(P.ratio , Clust_DP, Nzone, 
                                                  lwd=1.5, lty=2, col="grey25")}}
    #### Plot Algue Frac ####
    if(is.null(Malg_frac) == F){
      Mplot.alg <- Prep.data(Malg_frac, Seuil.Pour.alg, "Algue", Calc.stat = T, Only.keep = Keep.algue)
      
      if(Plot.NPP == T | Pol.influx.total == T){Clust_DP_plot = NULL}
      else(Clust_DP_plot = Clust_DP)
      
      if(Auto.Pollen.lab == T & is.null(Path.index) == F){
        InfoPol  <- read.csv(file = Path.index, sep=",",dec=".", header = T, stringsAsFactors = F, row.names = 1, check.names = FALSE)
        Algue.lab <- InfoPol$Label[match(gsub("\\.", " ", names(Mplot.alg[[1]])), InfoPol$Nom)]
      }
      else{Algue.lab <- NULL}
      Lab_algue_here <- Algue.lab 
      
      
      MAlgClean <- Mplot.alg[[1]]
      
      if(is.null(Keep.algue) == F & is.null(Order.algues) == T){
        Match.OK <- Keep.algue[which(Keep.algue %in% names(MAlgClean))]
        Add <- setdiff(names(MAlgClean), Match.OK)
        Save.order <- match(Add, names(MAlgClean))
        if(length(Add)>0){Match.OK <- c(Match.OK, Add)}
        Order.algues <- match(Match.OK, names(MAlgClean))
      }
      
      if(is.null(Order.algues) == F){
        MAlgClean <- MAlgClean[Order.algues]
        Lab_algue_here <- Lab_algue_here[Order.algues]
      }
      
      if(is.null(Malg_conc) == F & Al.Pol.influx == T){
        Axis.x.pos <- c(0, 1.4, 1.7)}
      else{Axis.x.pos <- NULL}
      
      p2 <- strat.plot(MAlgClean, yvar = as.double(Mplot.alg[[2]]), 
                       y.rev         = T, 
                       scale.percent = T,     # True pour des pourcentages
                       srt.xlabel    = 45,           # Rotation de 45 des noms de taxon
                       y.axis        = F,
                       xSpace        = 0.005,
                       x.names = Lab_algue_here,
                       x.pc.lab      = T,
                       x.axis = T,
                       cex.axis = 0.6,   # taille police légende en bas
                       #mgp = c(0, -0.7, -1), # position graduation en bas
                       mgp = Axis.x.pos, # horizon / vertical chiffre / vertical graduation
                       plot.poly     = T, 
                       plot.bar      = T,
                       plot.line     = T, 
                       lwd.bar       = 1,
                       col.poly      = "#a2dbe1ff",         # bleu clair
                       col.bar       = "#8c92a1ff",         # Fait apparaitre les traits continus
                       col.poly.line = "#8c92a1ff",         # gris fonce
                       clust = Clust_DP_plot,
                       clust.width=0.015, 
                       xLeft         = LAlgue,
                       xRight        = RAlgue,
                       add           = T
      )
      if(CONISS == TRUE & Nzone > 0){addClustZone(p2, Clust_DP, Nzone, 
                                                  lwd=1.5, lty=2, col="grey25")}
      
    }
    
    #### Plot Algue Conc. ####
    if(is.null(Malg_conc) == F & Al.Pol.influx == T){
      yo <- Prep.data(Malg_conc, Seuil.Conc.alg, "Algue", Only.keep = Keep.algue)
      spec = yo[[1]]
      Age = yo[[2]]
      couleur = yo[[3]]
      
      if(Plot.NPP == T | Pol.influx.total == T & is.null(Malg_frac) == F){Clust_DP_plot = NULL}
      else(Clust_DP_plot = Clust_DP)
      
      if(Auto.Pollen.lab == T & is.null(Path.index) == F){
        InfoPol  <- read.csv(file = Path.index, sep=",",dec=".", header = T, stringsAsFactors = F, row.names = 1, check.names = FALSE)
        Algue.lab <- InfoPol$Label[match(gsub("\\.", " ", names(Mplot.alg[[1]])), InfoPol$Nom)]
      }
      else{Algue.lab <- NULL}
      
      if(is.null(Malg_frac) == T){Lab_algue_here <- Algue.lab}
      else{Lab_algue_here = ""}
      
      if(is.null(Order.algues) == F){
        spec <- spec[Order.algues]
        if(is.null(Malg_conc) == T & Al.Pol.influx == T){ab_algue_here <- Lab_algue_here[Order.algues]}
      }  
      
      p2<- strat.plot(spec, yvar = as.double(Age), 
                      y.rev = TRUE, 
                      scale.percent = TRUE,     # True pour des pourcentages
                      srt.xlabel    = 45,       # Rotation de 45 des noms de taxon
                      y.axis        =FALSE,
                      xSpace        = 0.005,
                      x.pc.inc      = round(max(spec), digits = - 2)/4,
                      mgp = c(0, 0.5, 0.3),
                      plot.poly     = FALSE, 
                      x.names = Lab_algue_here,
                      plot.bar      = TRUE,
                      plot.line     = FALSE, 
                      lwd.bar       = 7,
                      col.poly      = "#a2dbe1ff",     # bleu clair
                      # col.bar       = couleur,         # Fait apparaitre les traits continus
                      col.bar       = "#004266",         # Fait apparaitre les traits continus
                      col.poly.line = "#8c92a1ff",     # gris fonce
                      exag=FALSE,                # Fait apparaitre la zone x10
                      col.exag="auto",           # Zone x10 couleur auto
                      clust = Clust_DP_plot,
                      clust.width=0.015, 
                      xLeft = LAlgue,
                      xRight = RAlgue,
                      add = TRUE
      )
      if(CONISS == TRUE & Nzone > 0){addClustZone(p2, Clust_DP, Nzone, 
                                                  lwd=1.5, lty=2, col="grey25")}}
    #### Plot NPP Conc. ####
    if(Plot.NPP == T){
      yo <- Prep.data(NPP_lake, Seuil.NPP, "NPP", Only.keep = Keep.NPP)
      spec = yo[[1]]
      #spec = spec/100
      Age = yo[[2]]
      couleur = yo[[3]]
      if(Spore.influx == T){spec["Spore Sum"] <- rowSums(spec)}
      if(Pol.influx.total == T | Diversite.plot == T){Clust_DP_plot = NULL}
      else(Clust_DP_plot = Clust_DP)
      if(Auto.Pollen.lab == T & is.null(Path.index) == F){
        InfoPol  <- read.csv(file = Path.index, sep=",",dec=".", header = T, stringsAsFactors = F, row.names = 1, check.names = FALSE)
        row.names(spec) <- gsub(" ",".", row.names(spec))
        InfoPol$Nom <- gsub(" ",".", InfoPol$Nom)
        NPP.lab <- InfoPol$Label[match(gsub("\\.", " ", names(spec)), gsub("\\.", " ", InfoPol$Nom))]
      }
      else{NPP.lab <- NULL}
      
      if(is.null(Keep.NPP) == F & is.null(Order.fungal) == T){
        Match.OK <- Keep.NPP[which(Keep.NPP %in% names(spec))]
        Add <- setdiff(names(spec), Match.OK)
        Save.order <- match(Add, names(spec))
        if(length(Add)>0){Match.OK <- c(Match.OK, Add)}
        Order.fungal <- match(Match.OK, names(spec))
      }
      
      if(is.null(Order.fungal) == F){
        spec <- spec[Order.fungal]
        NPP.lab <- NPP.lab[Order.fungal]
        couleur <- couleur[Order.fungal]
      }
      
      spec <- spec/10^3
      p3 <- strat.plot(spec, yvar = as.double(Age), 
                       y.rev         = TRUE, 
                       scale.percent = TRUE,     # True pour des pourcentages
                       srt.xlabel    = 45,           # Rotation de 45 des noms de taxon
                       y.axis        = FALSE,
                       xSpace        = 0.005, x.names = NPP.lab,
                       # x.pc.inc      = round(max(spec, na.rm = T), digits = - 2)/2,
                       x.pc.inc      = round(max(spec, na.rm = T), digits = - 2)/4,
                       mgp = c(0, 0.5, 0.3),
                       plot.poly     = FALSE, 
                       plot.bar      = TRUE,
                       plot.line     = FALSE, 
                       lwd.bar       = 7,
                       col.bar       = couleur,         # Fait apparaitre les traits continus
                       exag          = FALSE,                # Fait apparaitre la zone x10
                       col.exag      = "auto",           # Zone x10 couleur auto
                       clust         = Clust_DP_plot,
                       clust.width   = 0.015, 
                       xLeft         = LNPP,
                       xRight        = RNPP,
                       add           = TRUE
      )
      if(CONISS == TRUE & Nzone > 0){addClustZone(p3, Clust_DP, Nzone, 
                                                  lwd=1.5, lty=2, col="grey25")}}
    
    #### Plot Pollen inf. ####
    if(Pol.influx.total == T){
      if(is.null(Mpol_conc) == TRUE ){print("Calcul de l'influx pollinique impossible. Données manquantes.")}
      else{
        yo <- Prep.data(Mpol_conc, Seuil.Conc.alg, "Pollen")
        spec = yo[[1]]
        spec[["Pollen.tot.influx"]] <- rowSums(spec)/1000
        Age = yo[[2]]
        couleur = yo[[3]]
        
        if(Diversite.plot == T | Plot.algue == T | Plot.NPP == T ){Clust_DP_plot <-  NULL}
        else(Clust_DP_plot <- Clust_DP)
        P.inf <- strat.plot(spec$Pollen.tot.influx, yvar = as.double(Age), 
                            y.rev = TRUE, 
                            scale.percent = TRUE,     # True pour des pourcentages
                            srt.xlabel = 45,           # Rotation de 45 des noms de taxon
                            x.names = "Pollen Influx (#.10^3)",
                            y.axis=FALSE,
                            x.pc.inc      = round(max(spec$Pollen.tot.influx), digits = - 2)/2,
                            xSpace = 0.003,
                            mgp = c(0, 0.5, 0.3),
                            plot.poly     = FALSE, 
                            plot.bar      = TRUE,
                            plot.line     = FALSE, 
                            lwd.bar       = 7,
                            # col.poly      = "#a2dbe1ff",         # bleu clair
                            col.poly      = "#E6E6C6",          # gris clair
                            col.bar      = "#C67F05",          # gris clair
                            # col.bar       = couleur,         # Fait apparaitre les traits continus
                            col.poly.line = "#8c92a1ff",         # gris fonce
                            exag=FALSE,                # Fait apparaitre la zone x10
                            col.exag="auto",           # Zone x10 couleur auto
                            clust = Clust_DP_plot,
                            clust.width=0.015,
                            xLeft = LInf,
                            xRight = RInf,
                            add = TRUE
        )
        if(CONISS == TRUE & Nzone > 0){addClustZone(P.inf, Clust_DP, Nzone, 
                                                    lwd=1.5, lty=2, col="grey25")}}}
    
    #### Plot Diversity ####
    if(Diversite.plot == TRUE){
      MatDiv = spec[0]
      if(Pollen.diversity == T){
        yoPOL <- Prep.data(MP, Seuil.Pour, "Pollen")
        yoPOL <- data.frame(t(yoPOL[[5]]))
        yoPOL <- yoPOL[names(Age),]
        MatDiv["Pollen Diversity"] = yoPOL}
      
      if(Spore.diversity == T & Plot.NPP == T){
        YoNPP <- Prep.data(NPP_lake, Seuil.NPP)
        YoNPP <- data.frame(t(YoNPP[[5]]))
        YoNPP <- YoNPP[names(Age),]
        MatDiv["Spore Diversity"] = YoNPP}
      
      if(is.null(Ratio) == F){Clust_DP_plot = NULL}
      else{Clust_DP_plot = Clust_DP}
      
      pmat <- strat.plot(MatDiv, yvar = as.double(Age),
                         y.rev = TRUE, 
                         scale.percent = T,     # True pour des pourcentages
                         srt.xlabel = 45,           # Rotation de 45 des noms de taxon
                         #xlim = c(20,30),
                         xSpace = 0.003,
                         mgp = c(0, 0.5, 0.3),
                         y.axis        = F,
                         plot.poly     = T, 
                         plot.bar      = T,
                         plot.line     = T, 
                         lwd.bar       = 1,
                         col.poly      = "lightgrey",          # gris clair
                         col.bar       = "#8c92a1ff",         # Fait apparaitre les traits continus
                         col.poly.line = "#8c92a1ff",     # gris fonce
                         exag = F,                # Fait apparaitre la zone x10
                         col.exag = "auto",           # Zone x10 couleur auto
                         clust = Clust_DP_plot,
                         clust.width = 0.015, 
                         xLeft = LDiv,
                         xRight = RDiv,
                         add = T
      )
      if(CONISS == TRUE & Nzone > 0){addClustZone(pmat, Clust_DP, Nzone, 
                                                  lwd=1.5, lty=2, col="grey25")}
    }
    
    #### RETURN ####
    if(is.null(Save.plot) == F){dev.off()}
  }
  return(yo)
}

Plot.HYDE <- function(Extract, Models = c("baseline", "upper"), Limits = NULL, Log.scale = F, Facet_plots = F, 
                      Leg.pos = c(0.8,0.5),
                      Save.plot.RDS = NULL, X.lab.size = 8, Time.window = 500, Save.plot = NULL, W = 400, H = 400){
  #### Plot by country ####
  Extract.m <- melt(Extract, id = c("Age", names(Extract)[grep(paste(Models[[2]], "lower", sep = "|"), names(Extract))]))
  names(Extract.m)[names(Extract.m) == "variable"] <- "Mean"
  names(Extract.m)[names(Extract.m) == "value"] <- "Mean_value"
  
  Extract.m <- melt(Extract.m, id = c("Age",  "Mean", "Mean_value", names(Extract.m)[grep(Models[[2]], names(Extract.m))]))
  names(Extract.m)[names(Extract.m) == "variable"] <- "Low"
  names(Extract.m)[names(Extract.m) == "value"] <- "Low_value"
  
  Extract.m <- melt(Extract.m, id = c("Age", "Mean", "Mean_value", "Low", "Low_value"))
  Extract.m$color_mean <- gsub("_", "", gsub(paste(Models, collapse = "|"), "", Extract.m$Mean))
  Extract.m$color_up <- gsub("_", "", gsub(paste(c(Models, "lower"), collapse = "|"), "", Extract.m$variable))
  Extract.m$color_low <- gsub("_", "", gsub(paste(c(Models, "lower"), collapse = "|"), "", Extract.m$Low))
  
  Extract.m <- Extract.m[Extract.m$color_mean == Extract.m$color_low,]
  Extract.m <- Extract.m[Extract.m$color_mean == Extract.m$color_up,]
  Extract.m$facet <- ifelse(grepl("pop", Extract.m$color_mean), "Population counts", "Landscape uses (HYDE 3.2)")
  
  #### Param settings ####
  if(is.null(Limits) == T){Limits <- c(min(Extract.m$Age, na.rm = T), max(Extract.m$Age, na.rm = T))}
  
  if(Log.scale == T){
    To.log <- scale_x_log10()
    Xlabs <- xlab(expression(log(Surface)~(km^2)))}
  if(Log.scale == F){
    To.log <- NULL
    Xlabs <- xlab(expression(Surface~(km^2)))}
  
  My_color = c("cropland" = "#706432", "grazing" = "#7DA75E", "popc" = "darkred", "pasture" = "#1462BB")
  My_labs = c("cropland" = "Croplands", "grazing" = "Rangelands", "popc" = "Population counts")
  
  if(Facet_plots == T){My_facet <- facet_wrap(vars(facet), scale = "free_x")}
  if(Facet_plots == F){My_facet <- NULL}
  
  #### Plot ####
  Ptot <- ggplot(data = Extract.m, mapping = aes(y = Age, x = Mean_value, colour = color_mean))+ 
    geom_point()+ geom_line(orientation = "y")+ ylab("Age (yr cal BP)")+ Xlabs +
    My_facet +
    geom_colh(width = 25, color = NA, fill = NA, na.rm = T) +
    
    geom_ribbon(inherit.aes = F, mapping = aes(y = Age, xmin = value, xmax = Low_value, fill = color_mean), linetype = "dashed", alpha = .2)+
    scale_color_manual(values = My_color, labels = My_labs, name = "Human impact\n(HYDE 3.2)")+
    scale_fill_manual(values = My_color, labels = My_labs, name = "Human impact\n(HYDE 3.2)")+
    scale_y_reverse(limits = c(Limits[2], Limits[1]), breaks = seq(Limits[1], Limits[2], Time.window)) + To.log +
    theme(plot.background = element_blank(), panel.background = element_rect(fill = NA, colour = "grey30"),
          panel.grid = element_blank(), 
          axis.text.x = element_text(size = X.lab.size, angle = 45, hjust = 1, vjust = 1),
          axis.ticks.x = element_line(lineend = "butt", color = "grey70", linewidth = .1),
          
          legend.frame = element_blank(), legend.box = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), legend.position = Leg.pos)
  
  #### Export ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(Ptot, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.plot)}}
  
  if(is.null(Save.plot.RDS) == F){saveRDS(Ptot, Save.plot.RDS)}
  return(Ptot)
}

CWM.FT.plot <- function(MFT = NULL, MT = NULL, Mpol = NULL, Mbiom = NULL, Plot.x = NULL, Select.clim = NULL, Select.model = NULL, Multi.CWM = F, Zone.clim.box = T, Lab.age.rotation = 0,
                        Limites, Smooth.show = T, Biom.groups = T, Show.main.biomes = F, Eco.dot.size = 2.5, Eco.dot.alpha = .7, All.biomes = F, Taxa.rotation = 60, Panel.ecart = .1,
                        Select.DB = NULL, Select.trait = NULL, Strat.plot = T, Temp.zone, Select.interv = 1000, Smooth.param = 0.25, Pourc.DP = 25, Pourc.biom = 20, Zone.rotation = 0,
                        Manual.vlines = NULL, Nb.pol = 10, Plot.CWV = F, Taxa.hjust = 0.5, Trait.rotation = 90, Trait.hjust = 0.5, Pourc.lab.clim = 2,
                        Name.zone = NULL, Zone.clim = NULL, Xlims = NULL, Return.plot = F, Dot.size = 1, Dot.alpha = 0.1, DP.hist.size = 4, Select.pollen = NULL, Select.biom = NULL, X.title = NULL,
                        FT.select = NULL, R2.pos = NULL, Save.plot = NULL, Display.legend = "bottom", W = 1000, H = 1000, Facet_order = NULL){
  #### Settings ####
  if(missing(Temp.zone)){Temp.zone = rep("C", length(Zone.clim))}
  N.obs = 1; M.m1 <- NULL
  
  #### Temp zone ####
  if(any(unique(grepl("#", Temp.zone))) == T){
    print("Manual color scale for climate zone activated.")
    Rect.color.scale <- Temp.zone
    names(Rect.color.scale) <- Rect.color.scale
    if(any(unique(grepl("#", Temp.zone2))) == T){
      Rect.color.scale <- unique(c(Temp.zone, Temp.zone2))
      names(Rect.color.scale) <- Rect.color.scale
    }
    
  }
  else{
    Rect.color.scale <- c("W" = "#E76D51",
                          "C" = "#75AADB",
                          "G" = "grey40",
                          "D" = "#c67f05",
                          "Wt" = "#004266")
    Rect.color.scale <- Rect.color.scale[unique(Temp.zone)]
  }
  
  if(is.null(Zone.clim) == F){
    if(Zone.clim.box == T){Box.col <- "grey"}
    else{Box.col <- NA}
    
    yo <- data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                     xmax = Zone.clim[seq(2,length(Zone.clim), by=2)], 
                     Temp.col = Temp.zone)
    
    Zone.clim.to.plot <- geom_rect(data = yo, inherit.aes = F, na.rm = T,
                                   mapping = aes(ymin = -Inf, ymax= +Inf, xmin = xmin, xmax = xmax, fill = Temp.col),
                                   alpha=0.1, color = Box.col, linewidth = 0.3, linetype = 2)} 
  else{Zone.clim.to.plot <- NULL}
  
  if(is.null(Name.zone) == F){yo2 = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                                               xmax = Zone.clim[seq(2,length(Zone.clim), by=2)],
                                               Temp.col = Temp.zone,
                                               Title.zone = Name.zone,
                                               variable = "Height[CWM]"
  )}
  else{yo2 = data.frame(xmin = 0, xmax = 0, Temp.col = "", Title.zone = "", Categorie = "")}
  
  #### Clean Matrix function transfert ####
  if(is.null(MFT) == F){
    MFT <- data.frame(MFT)
    if(is.null(Select.clim) == T){
      print("Too much curves ! Please select *Select.clim*.")
      Select.clim <- c("MAAT", "MAP")
    }
    Match.clim <- unique (grep(paste(Select.clim,collapse="|"), names(MFT), value=TRUE))
    MFT <- MFT[Match.clim]
    
    if(is.null(Select.model) == F){
      Match.model <- unique (grep(paste(Select.model,collapse="|"), names(MFT), value=TRUE))
      MFT <- MFT[Match.model]
    }
    
    if(is.null(Select.DB) == F){
      Match.DB <- unique(grep(paste(Select.DB,collapse="|"), names(MFT), value=TRUE))
      MFT <- MFT[Match.DB]
    }
    MFT <- cbind(Site = row.names(MFT), MFT)
  }
  
  #### Clean Matrix trait ####
  if(is.null(MT) == F){
    if(Plot.CWV == T){MT.save <- MT}
    
    if(Multi.CWM == T){
      MT <- purrr::map(MT, 1)
      for(i in 1:length(MT)){MT[[i]]$Type <- names(MT)[i]}
      MT <- suppressMessages(purrr::reduce(MT, full_join))
    }
    else{
      if(class(MT) == "list"){MT <- MT[[1]]}
    }
    
    if(is.null(Plot.x) == T){
      Keep.x <- c("Age", "Top", "Bottom", "AgeBP2023", "Type")
      Plot.x <- names(MT)[names(MT) %in% Keep.x]
    }
    
    if(length(Plot.x)>0){
      if(Multi.CWM == F){Plot.x <- sort(Plot.x)[1]}
      if(Multi.CWM == T){Plot.x <- c("Age", "Type")}
    }
    if(is.null(Select.trait) == T){Select.trait <- names(MT)[grep("TRY_", names(MT))]}
    
    Select.trait <- c("Site", Plot.x, Select.trait)
    MT <- MT[Select.trait]
  }
  
  #### Clean Matrix trait CWV (variance) ####
  if(Plot.CWV == T){
    if(Multi.CWM == T){Mycondi <- "MCWV" %in% unlist(unique(purrr::map(MT.save, names)))}
    else{Mycondi <- "MCWV" %in% names(MT.save)}
    
    if(Mycondi == T){
      if(Multi.CWM == T){
        MV <- purrr::map(MT.save, 2)
        for(i in 1:length(MV)){MV[[i]]$Type <- names(MV)[i]}
        MV <- suppressMessages(purrr::reduce(MV, full_join))
      }
      else{
        if(class(MT.save) == "list"){MV <- MT.save[[2]]
        }
      }
      MV <- MV[Select.trait]
    }
    if(Mycondi == F){
      print("**** You want to plot the CWV, but this matrix is missing from MT list. ****")
      Plot.CWV = F
    }
  }
  
  #### Clean Matrix biomization ####
  if(is.null(Mbiom) == F){
    if(is.null(Plot.x) == T){
      Keep.x <- c("Age", "Top", "Bottom", "AgeBP2023")
      Plot.x <- names(Mbiom)[names(Mbiom) %in% Keep.x]
    }
    if(length(Plot.x)>0){Plot.x <- sort(Plot.x)[1]}
    
    if("BIOMPOL" %in% names(Mbiom) == T){Mbullet <- Mbiom[c(Plot.x, "BIOMPOL")]}
    if("BIOMPOL_" %in% names(Mbiom) == T){Mbullet <- Mbiom[c(Plot.x, "BIOMPOL_")]}
    
    if(is.null(Select.biom) == F){
      Select.biom <- c(Plot.x, Select.biom)
      Mbiom <- Mbiom[Select.biom]}
    else{
      if("BIOMPOL" %in% names(Mbiom) == T){Mbiom <- subset(Mbiom, select = -c(BIOMPOL))}
      if("BIOMPOL_" %in% names(Mbiom) == T){Mbiom <- subset(Mbiom, select = -c(BIOMPOL_))}
    }
    
    Mbiom <- melt(Mbiom, id = Plot.x)
    if(Biom.groups == F){Mbiom$groups <- "Biomization"}
    else{
      Mbiom$groups <- "Others"
      Mbiom$groups[grepl("CODE", Mbiom$variable)] <- "Deserts"
      Mbiom$groups[grepl("HODE", Mbiom$variable)] <- "Deserts"
      Mbiom$groups[grepl("ST", Mbiom$variable)] <- "Steppes"
      Mbiom$groups[grepl("X", Mbiom$variable)] <- "Woolands"
      Mbiom$groups[grepl("TEDE", Mbiom$variable)] <- "Forests~(decid.)"
      if(All.biomes == T){
        Mbiom$groups[grepl("TUND", Mbiom$variable)] <- "Tundra"
        Mbiom$groups[grepl("TAIG", Mbiom$variable)] <- "Forests~(Taiga/Conif.)"
        Mbiom$groups[grepl("COCO", Mbiom$variable)] <- "Forests~(Taiga/Conif.)"
        Mbiom$groups[grepl("CLDE", Mbiom$variable)] <- "Forests~(decid.)"
        Mbiom$groups[grepl("CLMX", Mbiom$variable)] <- "Forests~(mix.)"
        Mbiom$groups[grepl("COMX", Mbiom$variable)] <- "Forests~(mix.)"
      }
    }
  }
  #### Clean DP ####
  if(is.null(Mpol) == F){
    if(class(Mpol) == "list"){Mpol <- Mpol[[1]]}
    if(is.null(Plot.x) == T){
      Keep.x <- c("Age", "Top", "Bottom", "AgeBP2023")
      Plot.x <- names(Mpol)[names(Mpol) %in% Keep.x]
    }
    # if(any(names(Mpol) == "Site") == F){Mpol$Site <- row.names(Mpol)}
    if(length(Plot.x)>0){Plot.x <- sort(Plot.x)[1]}
    
    if(Plot.CWV == T){N.obs <- ncol(Mpol)}
    if(is.null(Select.pollen) == T){
      print("**** Only the 10 dominant taxa are displayed. You can manually choice them with 'Select.pollen', or change 10 to X with 'Nb.pol'.")
      Select.pollen <- Mpol[!names(Mpol) %in% Plot.x]
      Select.pollen <- names(colSums(Select.pollen)[order(colSums(Select.pollen), decreasing = T)])[1:Nb.pol]
    }
    
    Select.pollen <- c(Plot.x, Select.pollen)
    Mpol <- Mpol[Select.pollen]
    Mpol <- melt(Mpol, id = Plot.x)
    Mpol$Type <- "AP"
  }
  
  #### Select particular curves FT ####
  if(is.null(FT.select) == F & is.null(MFT) == F){
    FT.select <- c("Site", "Age", FT.select)
    MFT <- MFT[,which(names(MFT) %in% FT.select)]
  }
  
  #### Merge Matrixes ####
  if(is.null(MT) == F | is.null(MFT) == F){
    if(is.null(MT) == T){M <- MFT}
    if(is.null(MFT) == T){M <- MT}
    if(is.null(MFT) == F & is.null(MT) == F){M <- full_join(MT, MFT, by = intersect(names(MT), names(MFT)))}
    
    if("Site" %in% names(M) == T){M <- subset(M, select = -c(Site))}
    if(length(Plot.x) > 1){M <- M[order(M[[Plot.x[[1]]]]),]}
    else{M <- M[order(M[[Plot.x]]),]}
    
    M.m <- reshape2::melt(M, names(M)[which(!grepl("TRY_.", names(M)))])
    names(M.m)[c(ncol(M.m)-1,ncol(M.m))] <- c("Traits", "Traits_values")
    if(any(grepl("SEP",names(M.m)))){
      M.m <- reshape2::melt(M.m, names(M.m)[which(!grepl("SEP", names(M.m)))])
      names(M.m)[c(ncol(M.m)-1,ncol(M.m))] <- c("SEP", "SEP_values")
      M.m <- reshape2::melt(M.m, c(names(M.m)[1], "Traits", "Traits_values", "SEP", "SEP_values"))
    }
    else{
      M.m <- reshape2::melt(M.m, c(names(M.m)[1], "Traits", "Traits_values"))
    }
    
    if(is.null(MFT) == T){
      if(Multi.CWM == F){M.m$variable = NA; M.m$value = NA; M.m$Method = NA; M.m$Clim.param = NA}
    }
    else{
      M.m$Method <- sub("\\..*", "", M.m$variable)
      M.m$Clim.param <- sub(".*\\.", "", M.m$variable)
      
      M.m$variable <- gsub(paste(unique(M.m$Method), collapse = "|"), "", M.m$variable)
      M.m$variable <- gsub(paste(unique(M.m$Clim.param), collapse = "|"), "", M.m$variable)
      M.m$variable <- gsub("\\.", "", M.m$variable)
    }
    
    M.m$Traits <- gsub("TRY_", "", M.m$Traits)
    M.m$Traits[M.m$Traits == "LeafN"] = "N leaf"
    M.m$Traits[M.m$Traits == "LeafP"] = "P leaf"
    M.m$Traits[M.m$Traits == "LeafArea"] = "Leaf Area"
    M.m$Traits[M.m$Traits == "LeafThick"] = "Leaf Thickness"
    M.m$Traits[M.m$Traits == "Height"] = "Plant Height"
    M.m$Traits[M.m$Traits == "CNRatio"] = "C/N"
    # M.m$Traits[M.m$Traits == "SSD"] = "SSD[CWM]"
    # M.m$Traits[M.m$Traits == "SLA"] = "SLA[CWM]"
    M.m$Traits[M.m$Traits == "SeedMass"] = "Seed Mass"
  }
  else{M.m <- NULL}
  # print(names(M.m))
  # print(M.m)
  
  #### Add R2 ####
  if(Strat.plot == F & is.null(M.m) == F){
    if(R2.pos == "bottomleft"){
      R2.y = "bottom"
      R2.x = "left"}
    if(R2.pos == "bottomright"){
      R2.y = "bottom"
      R2.x = "right"}
    if(R2.pos == "none"){
      R2.y = "none"
      R2.x = "none"}
    
    R2.x <- 0.01
    # R2.x <- c(.95, .01, .95, .01, .01, .01)
    R2.y <- seq(from = 0.01, to = 0.5, length.out = 1+length(c(unique(M.m$variable), unique(M.m$Method))))
    
    # Add.r2 <- stat_poly_eq(label.y = R2.y, label.x = R2.x, color = "turquoise4", size = 3.5, small.r = F, na.rm = T,
    Add.r2 <- stat_poly_eq(label.y = R2.y, label.x = R2.x, size = 3, small.r = F, na.rm = T,
                           aes(label =  sprintf("%s*\", \"*%s" ,
                                                after_stat(rr.label),
                                                # after_stat(r.squared),
                                                after_stat(p.value.label)
                           )))
    
    if(length(unique(M.m$Clim.param)) == 1){
      Y.title <- unique(M.m$Clim.param)
      Y.title <- ylab(paste(Y.title, "(pollen-inferred)", sep = " "))
    }
    else{Y.title <- NULL}
  }
  
  #### Color settings ####
  Col.scale <- c("WASTDB" = "#e2a064ff", "WAST" = "#e2a064ff", "NMSDB" = "#F3A481",  "NAP" = "#c19d4dff", "AP" = "#72a72bff", "Pollen" = "grey40", "Full" = "grey40",
                 "ST" = "#91C4DD", "STDB" = "#91C4DD", "MDB" = "red", "ACADB" = "#963326",
                 "COSTDB" = "#c19d4dff", "COST" = "#c19d4dff", "TUSDB" = "#0094AF", "CAUCDB" = "#0F3361",
                 "MEDTEMP" = "#F3A481", "TEMPSCAND" = "#91C4DD", Traits = "#2c9740ff",
                 "EAPDB" = "#0F3361", "TAIGDB" = "#32156eff")
  
  My_colors <- scale_colour_manual(values = Col.scale, name = "Calibration databases")
  My_fill <- scale_fill_manual(values = Col.scale, name = "Calibration databases")
  
  # My_shape <- scale_shape(name = "Calibration databases")
  
  if(is.null(Xlims) == F){Xlims <- xlim(Xlims[1], Xlims[2])}
  if(is.null(Name.zone) == F & is.null(M.m) == F){yo2$variable = unique(M.m$Traits)[1]}
  if(is.null(Mbiom) == F | is.null(Mpol) == F | is.null(MT) == F){
    Zone.txt1 <- geom_text(data = yo2, inherit.aes = F, aes(x = (xmax+xmin)/2, y = 0, label = Title.zone, color = Temp.col), size = 3, angle = Zone.rotation, vjust = 0.5, hjust = 0.5, fontface = "bold")
  }
  else{Zone.txt1 <- NULL}
  
  if(is.null(Mbiom) == F){
    Couleur.Prentice <- c(
      "WAMX"="#185699FF",
      "TEDE"="#72a72bff",
      "XERO"="#cf0156ff",
      "COMX"="#c1b646ff",
      "HODE"="#dfd762ff",
      "WAST"="#e2a064ff",
      "CLMX"="#ee9b2fff",
      "PION"="#0e3056ff",
      "TAIG"="#32156eff",
      "TUND"="#caa6c2ff",
      "COCO"="#9d2e58ff",
      "COST"="#f3c768ff",
      "ANTH"="#6e1d2cff",
      "CODE"="#b0a625ff",
      "CLDE"="#b1d5f0ff",
      "No_data" = "grey80")
    
    Couleur.Prentice <- Couleur.Prentice[order(names(Couleur.Prentice))]
    Couleur.Prentice <- Couleur.Prentice[match(unique(Mbiom$variable), names(Couleur.Prentice))]
  }
  
  #### Plot RL ####
  if(Strat.plot == F){
    Plot.CWM <- ggplot(M.m, aes(y = value, x = Traits_values, colour = variable, shape = Method))+
      facet_wrap(vars(Traits), scales = "fixed") +
      geom_point(size = Dot.size, alpha = Dot.alpha)+
      geom_smooth(method = "lm", se = F, span = 1000, na.rm = T, linewidth = 1, formula = y ~ x)+
      Add.r2 + My_colors + Xlims + Y.title + xlab("CWM-traits (z-scores)")+
      guides(shape = guide_legend(override.aes = list(size = Dot.size+3))) +
      theme(
        plot.background = element_rect(fill = NA, colour = "grey30"),
        panel.spacing = unit(0.1, "lines"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside", panel.border = element_rect(NA, "black", linewidth = 1), strip.clip = "off",
        strip.text = element_text(size = 11),
        legend.justification = c("left"),               # left, top, right, bottom
        # plot.margin = unit(c(0,0,0,0), "lines")
      )
  }
  
  #### Plot Strat plot ####
  if(Strat.plot == T){
    #### Remelting ####
    if(is.null(M.m) == F){
      Keep.names.y <- names(M.m)[1]
      
      if(Multi.CWM == F){
        M1 <- M.m[c(1:3)]
        names(M1) <- c("xPlot", "variable", "value")
        M1$Method <- "Traits"
        M1$Clim.param <- "Traits"}
      else{
        M1 <- M.m
        names(M1) <- c("xPlot", "variable", "value", "Method", "Clim.param")
        M1$Method <- "Traits"
      }
      
      if(any(grepl("SEP",names(M.m)))){
        M2 <- M.m[c(1,9,7,8,6)]}
      else{
        if(Multi.CWM == F){
          M2 <- M.m[c(1,7,5,6,4)]
        }
        else{M2 <- M.m}
      }
      names(M2) <- c("xPlot", "variable", "value", "Method", "Clim.param")
      M.m <- rbind(M1,M2)
      
      # print(M.m)
      if(missing(Limites)){Limites = c(round(min(M.m$xPlot), digits = -2), max(M.m$xPlot))}
      
      M.m$variable[M.m$variable == "MAAT"] = "MAAT<br><span style='font-size:7pt; color:black'>(°C)</span>"
      M.m$variable[M.m$variable == "MTWAQ"] = "MTWAQ<br><span style='font-size:7pt; color:black'>(°C)</span>"
      M.m$variable[M.m$variable == "MTCO"] = "MTCO<br><span style='font-size:7pt; color:black'>(°C)</span>"
      M.m$variable[M.m$variable == "MTCOQ"] = "MTCOQ<br><span style='font-size:7pt; color:black'>(°C)</span>"
      M.m$variable[M.m$variable == "MAP"] = "MAP<br><span style='font-size:7pt; color:black'>(mm.yr<sup>-1</sup>)</span>"
      M.m$variable[M.m$variable == "MPCOQ"] = "MPCOQ<br><span style='font-size:7pt; color:black'>(mm.yr<sup>-1</sup>)</span>"
      M.m$variable[M.m$variable == "MPCO"] = "MPCO<br><span style='font-size:7pt; color:black'>(mm.yr<sup>-1</sup>)</span>"
      
      if(is.null(Facet_order) == T){Facet_order <- sort(unique(M.m$variable), decreasing = F)}
      
      M.m$variable <- factor(M.m$variable, levels = Facet_order, ordered = T)
    }
    else{
      if(is.null(Mpol) == F){Keep.names.y <- names(Mpol)[1]}
      if(is.null(Mbiom) == F){Keep.names.y <- names(Mbiom)[1]}
      if(is.null(MFT) == F){Keep.names.y <- names(MFT)[1]}
    }
    
    #### Settings ####
    if(Smooth.show == T){
      My_loess <- geom_smooth(method = "loess", se = F, fullrange = T, span = Smooth.param, na.rm = T, linewidth = .5, formula = 'y ~ x')
      My_lines <- NULL
    }
    else{My_loess <- NULL
    My_lines <- geom_line(linewidth = .1, na.rm = T)
    }
    
    if(is.null(Mbiom) == T & is.null(Mpol) == F & is.null(MT) == T & is.null(MFT) == T){Annot.Pol <- ""}
    if(is.null(Mbiom) == T & is.null(Mpol) == T & is.null(MT) == F & is.null(MFT) == T){Annot.CWM <- ""}
    if(is.null(Mbiom) == F & is.null(Mpol) == T & is.null(MT) == T & is.null(MFT) == T){Annot.Biom <- ""}
    if(is.null(Mbiom) == T & is.null(Mpol) == T & is.null(MT) == F & is.null(MFT) == F){Annot.CWM <- "(A)"; Annot.Clim <-"(B)"}
    if(is.null(Mbiom) == F & is.null(Mpol) == T & is.null(MT) == F & is.null(MFT) == T){Annot.Biom <-"(A)"; Annot.CWM <- "(B)"}
    if(is.null(Mbiom) == T & is.null(Mpol) == F & is.null(MT) == F & is.null(MFT) == T){Annot.Pol <- "(A)"; Annot.CWM <- "(B)"}
    if(is.null(Mbiom) == T & is.null(Mpol) == F & is.null(MT) == F & is.null(MFT) == F){Annot.Pol <- "(A)"; Annot.CWM <- "(B)"; Annot.Clim <-"(C)"}
    if(is.null(Mbiom) == F & is.null(Mpol) == T & is.null(MT) == F & is.null(MFT) == F){Annot.Biom <-"(A)"; Annot.CWM <- "(B)"; Annot.Clim <-"(C)"}
    if(is.null(Mbiom) == F & is.null(Mpol) == F & is.null(MT) == F & is.null(MFT) == T){Annot.Pol <- "(A)"; Annot.Biom <-"(B)"; Annot.CWM <- "(C)"}
    if(is.null(Mbiom) == F & is.null(Mpol) == F & is.null(MT) == F & is.null(MFT) == F){Annot.Pol <- "(A)"; Annot.Biom <-"(B)"; Annot.CWM <- "(C)"; Annot.Clim <- "(D)"}
    
    if(is.null(X.title) == T){
      if(Keep.names.y %in% c("Top", "Bottom", "Depth")){My.lab <- c("Depth (cm)")}
      if(Keep.names.y == "Age"){My.lab <- c("Time (cal. year BP)")}
    }
    else{My.lab <- X.title}
    
    My_theme <- theme(
      axis.title.y=element_text(size=12),
      axis.text.y = element_text(size = 8, colour = "grey30"),
      # axis.line.y = element_line(colour = "grey30"),
      axis.ticks.y = element_line(colour = "grey30"),
      legend.key.size = unit(5, "mm"),
      panel.spacing.x = unit(Panel.ecart, "lines"),
      legend.direction = "horizontal", plot.background = element_blank(),
      panel.spacing = unit(0.1, "lines"), legend.position = "bottom",
      panel.background = element_blank(), legend.background = element_blank(),
      strip.background = element_blank(), legend.spacing = unit(-0.6,"cm"),
      strip.placement = "left", strip.clip = "off",
      legend.justification = c("left"),
      plot.margin = unit(c(0,0,0,0), 'pt'))
    
    if(is.null(MFT) == F){My.axis.trait <- element_blank(); My.axis.trait2 <- element_blank(); My.axis.lab <- element_blank()}
    else{My.axis.trait <- NULL; My.axis.trait2 <- element_line(colour = "grey30", linewidth = .3); My.axis.lab <- element_text(angle = Lab.age.rotation, hjust = 1, vjust = 1)}
    
    if(is.null(MFT) == F | is.null(M.m) == F){My.axis.biom <- element_blank(); My.axis.biom2 <- element_blank()}
    else{My.axis.biom <- NULL; My.axis.biom2 <- element_line(colour = "grey30", linewidth = .3)}
    
    #### Plot CWM ####
    if(is.null(M.m) == F){
      #### Clean ####
      M.m1 <- M.m[M.m$Method == "Traits",]
      M.m2 <- M.m[M.m$Method != "Traits",]
      M.m1 <- M.m1[!duplicated(M.m1),]
      M.m2 <- M.m2[!duplicated(M.m2),]
      M.m1$variable <- factor(as.character(M.m1$variable))
      M.m1 <- M.m1[!is.na(M.m1$variable),]
      M.m1 <- M.m1[M.m1$xPlot <= Limites[2],]
      M.m1 <- M.m1[M.m1$xPlot >= Limites[1],]
      
      #### Add ribbon from CWV ####
      if(Plot.CWV == T){
        if("Site" %in% names(MV) == T){MV <- subset(MV, select = -c(Site))}
        
        if(length(Plot.x) > 1){MV <- MV[order(MV[[Plot.x[[1]]]]),]}
        else{MV <- MV[order(MV[[Plot.x]]),]}
        
        MV.m <- reshape2::melt(MV, names(MV)[which(!grepl("TRY_.", names(MV)))])
        
        MV.m$variable <- gsub("TRY_", "", MV.m$variable)
        MV.m$variable[MV.m$variable == "LeafN"] = "N leaf"
        MV.m$variable[MV.m$variable == "LeafP"] = "P leaf"
        MV.m$variable[MV.m$variable == "CNRatio"] = "C/N"
        MV.m$variable[MV.m$variable == "LeafArea"] = "Leaf Area"
        MV.m$variable[MV.m$variable == "LeafThick"] = "Leaf Thickness"
        MV.m$variable[MV.m$variable == "Height"] = "Plant Height"
        # MV.m$variable[MV.m$variable == "SSD"] = "SSD[CWM]"
        # MV.m$variable[MV.m$variable == "SLA"] = "SLA[CWM]"
        MV.m$variable[MV.m$variable == "SeedMass"] = "Seed Mass"
        names(MV.m)[names(MV.m) == "value"] <- "variance"
        names(MV.m)[names(MV.m) == "Type"] <- "Clim.param"
        names(MV.m)[1] <- "xPlot"
        
        if(Multi.CWM == F){MV.m$Clim.param <- "Traits"}
        
        M.m1 <- full_join(M.m1, MV.m, by = c("xPlot", "Clim.param", "variable"))
        
        M.m1$MinI <- M.m1$value - 1.96*sqrt(M.m1$variance)/sqrt(N.obs)
        M.m1$MaxI <- M.m1$value + 1.96*sqrt(M.m1$variance)/sqrt(N.obs)
        
        Envellope <- geom_ribbon(mapping = aes(ymin = MinI, ymax = MaxI,  fill = Clim.param), colour = NA, linewidth = 0, alpha = .3)
        My_fill <- My_fill
        New.scale <- new_scale_fill()
      }
      else{Envellope <- NULL; My_fill <- NULL; New.scale <- NULL}
      
      #### Final Plot ####
      Plot.CWM <- ggplot(M.m1, aes(y = value, x = xPlot, colour = Clim.param, shape = Method))+
        scale_fill_manual(values = Rect.color.scale, guide = "none")+
        Zone.clim.to.plot + 
        New.scale +
        My_fill + #My_shape +
        Envellope +
        facet_wrap(vars(variable), scales = "free_y", ncol = 1, drop = F, strip.position = "left") + xlab(My.lab)+
        # facet_wrap(vars(variable), scales = "free_y", ncol = 1, drop = F, strip.position = "left", labeller = label_parsed) + xlab(My.lab)+
        geom_point(size = Dot.size, alpha = Dot.alpha)+ labs(y = paste(Annot.CWM, "CWM traits"))+
        My_lines + My_loess + My_colors + 
        scale_x_continuous(breaks = c(Limites[1], round(seq(0, Limites[2], by = Select.interv))))+
        scale_y_continuous(breaks = scales::breaks_extended(n = 4))+
        geom_rangeframe(sides = "l", color = "grey30", size = .3)+
        coord_cartesian(clip="off") +
        geom_vline(xintercept = Manual.vlines, col = "grey30", lty = 2, alpha = 0.7)+
        guides(shape = "none", colour = guide_legend(nrow = 1, override.aes=list(fill=NA, size = Dot.size*2, linewidth = 0.7)))+
        theme(
          strip.background = element_blank(),
          panel.spacing.y = unit(Panel.ecart, "lines"),
          strip.text.y.left = element_text(angle = Trait.rotation, hjust = Trait.hjust),
          axis.title.x = My.axis.trait,
          axis.ticks.x = My.axis.trait,
          axis.line.x = My.axis.trait2,
          legend.background = element_blank(),
          axis.text.x = My.axis.lab
        )}
    else{Plot.CWM <- NULL}
    
    #### Plot FT ####
    if(is.null(MFT) == F){
      # M.m1$variable <- as.factor(as.character(M.m1$variable))
      M.m2$variable <- as.factor(as.character(M.m2$variable))
      M.m2 <- M.m2[M.m2$xPlot <= Limites[2],]
      M.m2 <- M.m2[M.m2$xPlot >= Limites[1],]
      
      Plot.Clim <- ggplot(M.m2, aes(y = value, x = xPlot, colour = Clim.param, shape = Method))+
        scale_fill_manual(values = Rect.color.scale, guide = "none")+
        Zone.clim.to.plot +
        facet_wrap(vars(variable), scales = "free_y", ncol = 1, drop = F, strip.position = "left") + xlab(My.lab)+
        geom_point(size = Dot.size, alpha = Dot.alpha)+
        # geom_point()+ 
        labs(y = paste(Annot.Clim, "Climate reconstructions"))+
        My_lines + My_loess + My_colors + #My_shape +
        scale_x_continuous(breaks = c(Limites[1], round(seq(0, Limites[2], by = Select.interv))))+
        scale_y_continuous(expand = c(0.05,0.05))+
        # new_scale_color()+
        # scale_color_manual(values = Rect.color.scale, guide = "none", name = NULL, labels = NULL, breaks = NULL, na.translate = FALSE)+
        # Zone.txt1 +
        geom_vline(xintercept = Manual.vlines, col = "grey30", lty = 2, alpha = 0.7)+
        
        guides(shape = guide_legend(nrow = 1, override.aes = list(size = Dot.size*2)), 
               colour = guide_legend(nrow = 1, override.aes=list(fill=NA, size = Dot.size*2, linewidth = 0.7)))+
        theme(
          axis.title.x=element_text(size=12, colour = "grey20"),
          axis.line.x = element_line(colour = "grey30"), 
          legend.background = element_blank(),
          axis.ticks.x = element_line(colour = "grey30"), strip.text = element_markdown(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "grey30"))
      
      
      Plot.CWM <- Plot.CWM / Plot.Clim & My_theme
    }
    else{
      Plot.CWM <- Plot.CWM & My_theme
      
    }
    
    #### Plot Biomization ####
    if(is.null(Mbiom) == F){
      library(tidypaleo)
      names(Mbiom)[1] <- "xPlot"
      if(missing(Limites)){Limites = c(round(min(Mbiom$xPlot, na.rm = T), digits = -2), max(Mbiom$xPlot, na.rm = T))}
      
      Mbiom <- Mbiom[Mbiom$xPlot <= Limites[2],]
      Mbiom <- Mbiom[Mbiom$xPlot >= Limites[1],]
      
      if(nlevels(Mbiom$variable) >= 6){Nrow.biomes = 2}
      if(nlevels(Mbiom$variable) > 10){Nrow.biomes = 3}
      else{Nrow.biomes = 1}
      
      Plot.biom <- ggplot(Mbiom, aes(y = value, x = xPlot, color = variable)) +
        scale_fill_manual(values = Rect.color.scale, guide = "none")+
        Zone.clim.to.plot + 
        facet_wrap(vars(groups), scales = "free_y", ncol = 1, drop = F, strip.position = "left", labeller = label_parsed) + xlab(My.lab)+
        geom_point(size = Dot.size, alpha = Dot.alpha, shape = 15)+
        My_loess + My_lines +
        labs(y = paste(Annot.Biom, "Biome-scores"))+
        scale_x_continuous(breaks = c(Limites[1], round(seq(0, Limites[2], by = Select.interv))))+
        
        scale_y_continuous(breaks = scales::breaks_extended(n = 4))+
        geom_rangeframe(sides = "l", color = "grey30", size = .3)+
        coord_cartesian(clip="off") +
        scale_color_manual(values = Couleur.Prentice, label = names(Couleur.Prentice), name = "Biomes")+
        geom_vline(xintercept = Manual.vlines, col = "grey30", lty = 2, alpha = 0.7)+
        
        guides(shape = "none", colour = guide_legend(nrow = Nrow.biomes, override.aes=list(fill=NA, size = Dot.size*2, linewidth = 0.7)))+
        theme(
          axis.title.x = My.axis.biom,
          axis.ticks.x = My.axis.biom,
          axis.line.x = My.axis.biom2,
          axis.text.x = My.axis.biom,
          # axis.line.y = element_line(colour = "grey30"),
          legend.key.size = unit(5, "mm"),
          legend.direction = "horizontal", legend.key.spacing.y = unit(0,"cm"),
          panel.spacing = unit(Panel.ecart, "lines"),
          legend.spacing = unit(-0.6,"cm"), 
          panel.background = element_blank(), legend.background = element_blank(),
          strip.background = element_blank(), plot.background = element_blank(),
          strip.placement = "left", strip.clip = "off",
          plot.margin = unit(c(0,0,0,0), 'pt')
        )
      
      if(is.null(M.m) == F){
        Plot.CWM <- Plot.biom/Plot.CWM/guide_area() + plot_layout(guides = "collect", heights = c(Pourc.biom, 100-Pourc.biom, 10)) & 
          theme(legend.spacing = unit(-0.4,"cm"), plot.margin = unit(c(0,0,0,0), "mm"))}
      if(is.null(M.m) == T){Plot.CWM <- Plot.biom&My_theme}
    }
    
    #### Plot DP ####
    if(is.null(Mpol) == F){
      library(tidypaleo)
      
      names(Mpol)[1] <- "xPlot"
      if(missing(Limites)){Limites = c(round(min(Mpol$xPlot, na.rm = T), digits = -2), max(Mpol$xPlot, na.rm = T))}
      
      Mpol <- Mpol[Mpol$xPlot <= Limites[2],]
      Mpol <- Mpol[Mpol$xPlot >= Limites[1],]
      
      Plot.DP <- ggplot(Mpol, aes(y = value*100, x = xPlot)) +
        scale_fill_manual(values = Rect.color.scale, guide = "none")+
        Zone.clim.to.plot + 
        new_scale_fill()+
        scale_fill_manual(values = c("AP" = "#c19d4dff"), guide = "none")+
        geom_ribbon(aes(ymin = 0, ymax = value*100, fill = Type), position = "identity", alpha = .3) +
        geom_hline(yintercept = 0, color = "grey", size = 0.7)+
        geom_line(color = "grey30", orientation = "x")+
        geom_col(width = DP.hist.size, position = "dodgev", na.rm = T) +
        facet_grid(variable ~ ., space = "free_y", scales = "free", switch = "y")+
        labs(y = paste(Annot.Pol, "Pollen FA (%)"))+
        scale_y_continuous(breaks = seq(10,100,10), expand = c(0,0))+
        geom_vline(xintercept = Manual.vlines, col = "grey30", lty = 2, alpha = 0.7)+
        
        theme(
          # axis.title.y = element_text(size=12),
          axis.title.x = element_blank(),
          # strip.text.y = element_text(size = 8, angle = 145),
          strip.text.y.left = element_markdown(size = 8, angle = Taxa.rotation, hjust = Taxa.hjust),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8, colour = "grey30"),
          axis.line.y = element_line(colour = "grey30"),
          axis.line.x = element_blank(),
          axis.ticks.y = element_line(colour = "grey30"),
          axis.ticks.x = element_blank(),
          legend.position = "none", 
          legend.key.size = unit(8, "mm"),
          legend.direction = "horizontal", #
          panel.spacing = unit(Panel.ecart, "lines"), plot.background = element_blank(),
          panel.background = element_blank(), legend.background = element_blank(),
          strip.background = element_blank(),
          strip.placement = "left", strip.clip = "off",
          legend.justification = c("left"),
          plot.margin = unit(c(0,0,0,0), 'pt')
        )
      
      if(is.null(Plot.CWM) == T){Plot.CWM <- Plot.DP/Plot.CWM}
      else{Plot.CWM <- Plot.DP/Plot.CWM + plot_layout(heights = c(Pourc.DP, 100-Pourc.DP)) }
      
    }
    
    #### Plot biomes as bullet ####
    if(Show.main.biomes == T & is.null(Mbiom) == F){
      names(Mbullet)[names(Mbullet) != Plot.x] <- "Biome"
      names(Mbullet)[names(Mbullet) == Plot.x] <- "variable"
      
      Mbullet <- Mbullet[Mbullet$variable <= Limites[2],]
      Mbullet <- Mbullet[Mbullet$variable >= Limites[1],]
      
      Peco  <-  ggplot(Mbullet, mapping = aes(y = 1, x = variable, colour = Biome))+
        geom_point(size = Eco.dot.size, shape = 15, na.rm = T, alpha = Eco.dot.alpha)+
        scale_color_manual(values = Couleur.Prentice, label = names(Couleur.Prentice), name = "Biomes")+
        coord_cartesian(clip = 'off') +
        theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
              panel.background=element_blank(), panel.border = element_blank(),
              plot.background = element_blank(), panel.spacing = unit(Panel.ecart, 'lines'),
              plot.margin = unit(c(0,0,0,0), "lines"),
              legend.position = "none", strip.clip = "off",
              panel.grid = element_blank())
      
      Plot.CWM <- Peco/Plot.CWM +  plot_layout(heights = c(2, 98))
    }
    
    #### Plot Temp zones ####
    if(is.null(Name.zone) == F & is.null(M.m1) == F){
      M.annot <- M.m1[c(which.min(M.m1$xPlot), which.max(M.m1$xPlot)),]
      M.annot$value <- 0
      
      Plot.annot <- ggplot(M.annot, aes(y = value, x = xPlot))+
        geom_point(colour = "white", alpha = 0.01)+
        scale_fill_manual(values = Rect.color.scale, guide = "none")+
        labs(y = paste(Annot.CWM, "CWM traits"))+ 
        Zone.txt1 +
        coord_cartesian(clip="off") +
        scale_color_manual(values = Rect.color.scale, guide = "none", name = NULL, labels = NULL, breaks = NULL, na.translate = FALSE)+
        theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank(),
              strip.background = element_blank(), panel.grid = element_blank(), legend.background = element_blank(), 
              axis.title = element_blank(), panel.spacing = unit(Panel.ecart, 'lines'),
              axis.ticks = element_blank(), strip.clip = "off",
              plot.margin = unit(c(0,0,0,0), 'pt'),
              axis.text = element_blank()
        )
      
      Plot.CWM <- Plot.annot/Plot.CWM +  plot_layout(heights = c(Pourc.lab.clim, 100-Pourc.lab.clim)) & theme(panel.grid = element_blank(), legend.background = element_blank())
    }
  }
  
  
  #### Save plot and export ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(Plot.CWM, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.plot)}}
  
  if(Return.plot == F){return(M.m)}
  else{return(Plot.CWM)}
  
}

Compar.Surf.Calib <- function(MGDGT, MClim, Meco, Select.clim, Clim.display, Calib.to.keep = NULL, Leg.row, H.arrow, Box.order.mean = F, Add.frame = F, Break.line.annot = T,
                              Order.by, Bold_dots = NULL, Barres, Lab.angle = 45, Lab.size = 11, Bold.width = 3, Arrow.width = 0.5, Lab.area, Units, Y.limits, Save.plot, H, W, Pt.size){
  #### Settings ####
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Clim.display)){Clim.display = "Normal"}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Pt.size)){Pt.size = 1}
  if(missing(Select.clim)){Select.clim = "MAAT"}
  if(missing(Units)){Units = "(°C)"}
  if(missing(Leg.row)){Leg.row = NULL}
  if(missing(Barres)){Barres = NULL}
  if(missing(Lab.area)){Lab.area = NULL}
  if(missing(Meco)){Meco = NULL}
  if(missing(Order.by)){Order.by = NULL}
  if(missing(H.arrow)){H.arrow = F}
  
  #### Select model ####
  if(is.null(Calib.to.keep) == F){Clim.modelled <- MGDGT[names(MGDGT)%in% Calib.to.keep]}
  Clim.modelled <- Clim.modelled[,grep(paste("^", Select.clim, ".", sep = ""), colnames(Clim.modelled))]
  
  Clim.real <- MClim[[Select.clim]]
  Order.site <- factor(row.names(MGDGT))
  
  if(Clim.display == "Normal"){
    T.ordonnee = paste(Select.clim, Units, sep = " ")
    Mat.graph <- cbind(Clim.modelled, Clim.real)}
  
  if(Clim.display == "Anomaly"){
    T.ordonnee = bquote(Delta ~.(Select.clim) ~ .(Units))
    Mat.graph <- Clim.real - Clim.modelled}
  
  if(Clim.display == "Error"){
    T.ordonnee = bquote(Delta~.(Select.clim) ~ "(% error)")
    Mat.graph <- (Clim.real - Clim.modelled)/Clim.real}
  
  if(Clim.display == "Pourcentage.error"){
    T.ordonnee = bquote(Delta~.(Select.clim) ~ "(% error)")
    Mat.graph <- (Clim.real - Clim.modelled)/Clim.real*100}
  
  if(missing(Y.limits)){Y.limits = c(min(Mat.graph), max(Mat.graph))}
  
  if(Box.order.mean == T){
    New_order <- order(abs(colMeans(Mat.graph)))
    Order_mean <- names(Mat.graph)[New_order]
    My_arrow_lab <- "Decreasing mean"
    
  }
  else{My_arrow_lab <- "Increasing param nb"}
  
  #### Melt + graphical settings ####
  Mat.graph <- cbind(Mat.graph, Sites = rownames(MGDGT))
  
  if(is.null(Meco) == F){
    names(Meco)[1] <- "Clustering"
    Meco$Sites <- row.names(Meco)
    Merror_graph <- dplyr::left_join(Mat.graph, Meco)
    Merror_graph <- melt(Merror_graph, id =c('Sites', 'Clustering'))
  }
  else{
    Merror_graph <- melt(Mat.graph, id ='Sites')}
  
  if(is.null(Order.by) == F){
    New.order <- row.names(Order.by)[order(Order.by[[1]])]
    Merror_graph$Sites <- factor(Merror_graph$Sites, levels = New.order, ordered = T)
  }
  
  Model.lab <- paste(sub("_", "[", colnames(Mat.graph)[-ncol(Mat.graph)]), "]", sep = "")
  Model.lab <- gsub("_", "~", Model.lab)
  Model.lab <- gsub("5", "5*", Model.lab)
  
  if(is.null(Barres) == F){
    Barre1 = Barres[1]
    if(is.na(Barres[2])){Barre2 <- nrow(MGDGT)+1}
    else{Barre2 = Barres[2]}
    
    Barre3 = Barres[3]
  }
  else{Barre1 = NULL; Barre2 = NULL; Barre3 = NULL}
  
  Label1 = Lab.area[1]
  Label2 = Lab.area[2]
  Label3 = Lab.area[3]
  Label4 = Lab.area[4]
  
  if(Select.clim %in% c("MAAT", "MAF")){Lab.up <- "Over\nwarm"; Lab.down <- "Over\ncold"}
  if(Select.clim %in% c("MAP", "MPCOQ", "AI")){Lab.up <- "Over\ndry"; Lab.down <- "Over\nwet"}
  
  if(Break.line.annot == F){Lab.up <- gsub("\n", " ", Lab.up); Lab.down <- gsub("\n", " ", Lab.down)}
  if(Add.frame == T){
    My_frame <- element_rect(fill = NA, color = "black")
    My_axes <- element_blank()}
  else{
    My_frame <- element_blank()
    My_axes <- element_line(colour = "black")}
  
  #### Bold dots ####
  if(is.null(Bold_dots) == F){
    Mbold <- Merror_graph[Merror_graph$variable == Bold_dots,]
    Add_bold <- geom_dotplot(data = Mbold, binaxis = "y", binwidth = Pt.size, alpha = 0.7, na.rm = T, fill = NA, color = "black", stroke = Bold.width)
  }
  else{Add_bold <- NULL}
  
  #### Arrows ####
  if(H.arrow == T){
    Horiz.arrow <- geom_segment(x = 0.3, y = Y.limits[2], xend = length(unique(Merror_graph$variable))-2, yend = Y.limits[2],
                                arrow = arrow(length = unit(0.02, "npc"), ends = "last"),
                                arrow.fill = "grey30", colour = "grey30", size = Arrow.width)
    Horiz.annot <- annotate("text", x = 0.3, y = Y.limits[2]-Y.limits[2]*.1, hjust = 0, label = My_arrow_lab, size = 4, colour = "grey30")
  }
  else{
    Horiz.arrow <- NULL
    Horiz.annot <- NULL}
  
  #### Order boxes by mean values ####
  if(Box.order.mean == T){
    Merror_graph$variable <- factor(Merror_graph$variable, levels = rev(Order_mean), ordered = F)
    Model.lab <- Model.lab[rev(New_order)]
  }
  else{Merror_graph$variable <- factor(Merror_graph$variable, levels = levels(Merror_graph$variable), ordered = F)}
  
  #### Plot Clim error matrice ####
  p1 <- ggplot(data = Merror_graph , aes(x = Sites, y = value, fill = variable)) +
    geom_abline(slope = 0, intercept = 0, lty="dashed")+
    geom_dotplot(binaxis = "y", binwidth = Pt.size, alpha = 0.7, na.rm = T, color = NA) +
    Add_bold + 
    xlab("Surface samples") +                                                                                                 # Légende abscisse
    ylab(T.ordonnee)+                                                                                               # Légende ordonnée
    # scale_x_discrete(limits = Order.site) +
    scale_shape_manual(values = c(1,1,1,1,1,17,21,15))+
    scale_fill_discrete(labels = parse(text = Model.lab))+
    geom_vline(xintercept = Barre1, lty="dotted")+
    geom_vline(xintercept = Barre2, lty="dotted")+
    geom_vline(xintercept = Barre3, lty="dotted")+
    ggplot2::annotate("text", x = Barre1/2, y = Y.limits[2], label = Label1, size = 5)+
    ggplot2::annotate("text", x = (Barre1 + Barre2)/2, y = Y.limits[2], label = Label2, size = 5)+
    ggplot2::annotate("text", x = (Barre2 + Barre3)/2, y = Y.limits[2], label = Label3, size = 5)+
    ggplot2::annotate("text", x = (Barre3 + nrow(MGDGT))/2, y = Y.limits[2], label = Label4, size = 5)+
    ylim(Y.limits[1],Y.limits[2])+
    guides(fill = guide_legend(override.aes = list(size = Pt.size*20, alpha = 1), nrow = Leg.row))+
    
    theme(
      axis.title = element_text(margin = ggplot2::margin(t = -2, r = 0, b = 0, l = 0), size = 15),
      axis.text.x = element_text(angle = Lab.angle, hjust = 1, size = Lab.size),
      axis.text.y = element_text(hjust = 1, size = 14),
      plot.background = element_blank(), panel.grid = element_blank(),
      panel.background = My_frame, panel.border = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",                   # permet de mettre le carre des legendes en bas
      legend.direction = "horizontal",
      legend.text.align = 0,
      legend.text = element_text(size = 14),
      legend.key = element_blank(),
      axis.line = My_axes
    )
  
  #### Plot boxplot #####
  if(is.null(Meco) == F){
    p2 <- ggplot(Merror_graph, aes(x = variable, y = value, fill = Clustering, color = variable))+
      guides(color = "none")
  }
  else{
    p2 <- ggplot(Merror_graph, aes(x = variable, y = value, fill = as.factor(variable)))+
      guides(fill = "none")+
      geom_jitter(width = 0.1, shape = 1, alpha = 0.7, na.rm = T)
  }
  
  #### Plot features ####
  p2 <- p2 + 
    geom_boxplot(outlier.colour = "red", outlier.shape = NA, alpha = 0.6, na.rm = T, notch = F) +
    xlab(paste(Select.clim, "calibrations", sep = " ")) +                                                                                                 # Légende abscisse
    scale_x_discrete(labels = parse(text = Model.lab))+
    geom_abline(slope = 0, intercept = 0, lty="dashed")+
    
    ylim(Y.limits[1],Y.limits[2])+
    #### Annotations ####
  geom_segment(x = length(unique(Merror_graph$variable))+0.5, y = (Y.limits[1] - Y.limits[1]*0.15), xend = length(unique(Merror_graph$variable))+0.5, yend = (0+Y.limits[1]*0.15),
               arrow = arrow(length = unit(0.02, "npc"), ends = "first"),
               arrow.fill = "#A30A2A", colour = "#A30A2A", size = Arrow.width)+
    geom_segment(x = length(unique(Merror_graph$variable))+0.5, y = (Y.limits[2] - Y.limits[2]*0.15), xend = length(unique(Merror_graph$variable))+0.5, yend = (0+Y.limits[2]*0.15),
                 arrow = arrow(length = unit(0.02, "npc"), ends = "first"),
                 arrow.fill = "#3B6CAF", colour = "#3B6CAF", size = Arrow.width)+
    annotate("text", x = length(unique(Merror_graph$variable))+0.3, y = Y.limits[2]*.95, label = Lab.down, size = 3, colour = "#3B6CAF", hjust = 1)+
    annotate("text", x = length(unique(Merror_graph$variable))+0.3, y = Y.limits[1]*.95, label = Lab.up, size = 3, colour = "#A30A2A", hjust = 1)+
    Horiz.arrow +
    Horiz.annot + 
    
    #### Theme ####
  theme(
    axis.title.x = element_text(margin = ggplot2::margin(-2, 0, 0, 0), size = 15),
    axis.text.x  = element_text(angle = Lab.angle, hjust = 1, size = Lab.size),
    plot.background = element_blank(), panel.grid = element_blank(),
    panel.background = My_frame, panel.border = element_blank(),
    axis.line.x = My_axes, 
    axis.line.y  = element_blank() ,           # fait apparaitre seulement l'axe des x
    axis.text.y  = element_blank(),            # cache les valeurs de l'axe y
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()#,           # cache les graduations sur l'axe
  )
  
  #### Return ####
  pf <- p1 + p2 + plot_layout(ncol = 2, nrow = 1, widths = c(2/3, 1/3)) 
  pf <- pf + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 18, vjust = 1.3, hjust = -0.3))
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){
      ggsave(pf, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.plot)}}
  return(pf)
}


MGDGT.change.name <- function(List.GDGT = NULL, Old.names = NULL, New.names = NULL) {
  if (length(New.names) != length(Old.names)) {stop("Old and new names vectors must have the same length and order.")}
  renamed_list <- lapply(List.GDGT, function(df) {for (i in seq_along(Old.names)) {if (Old.names[i] %in% names(df)) {names(df)[names(df) == Old.names[i]] <- New.names[i]}}
    return(df)})
  return(renamed_list)}

Plot.GDGT.clim <- function(MGDGT, MMin, MMax, Pclim, Select.interv, Surf.val, Anomaly, Title.age, CRUTS, Manual.vlines = NULL, Manual.color.scale = NULL,
                           Facet.T, Facet.scale.egal, Cores.lab, Clim.lab, Select.model, Smooth.SD, Display.legends = NULL, Panel.box = F, RMSE.barres = NULL,
                           Label.group, Group.name, Age.select, Remove.models, Manual.y.val, Zoom.box, Show.repels = T, Repel.repoussage = 10,
                           Show.x.axis = T, Shape.group = NULL, Keep.models = NULL, Dot.size = 2, Lab.clim.angle = 0, Surf.val.max.age = NULL,
                           Pollen.plot.merge, Mono.core, Multi.clim, Pourc.clim.param, Hiatus.age, Temp.col.alpha = 0.1, Color.by.proxy = F,
                           Keep.other.grp, Annotation, Facet.size.egal, Smooth.param, Extract.fitting, Repel.x = NULL, Zone.clim.dashed = T,
                           Detrending.fitting = F,
                           Save.plot, W, H, Save.Rds, Limites, Add.lim.space, Zone.clim, Name.zone, Temp.zone){
  #### Libraries ####
  if(missing(MGDGT)){warning("Import the GDGT matrix for plotting.")}
  if(missing(Pclim)){warning("Select and sort the climate parameter(s) to plot.")}
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(reshape)      # permet d'utiliser la fonction melt utile pour afficher les graphs
  library(cowplot)      # permet aligner les graphs  
  library(purrr) # transpose les matrices
  library(ggrepel) # nom des lignes à côté
  library(ggnewscale) # function new_scale_color
  library(RColorBrewer)
  library(stringr) # str_count
  library(patchwork)
  library(zoo) # use function na.locf
  library(ggthemes) # geom_rangeframe
  
  #### Init Val ####
  if(missing(Hiatus.age)){Hiatus.age = NULL}
  if(missing(MMin)){MMin = NULL}
  if(missing(MMax)){MMax = NULL}
  if(missing(Surf.val)){
    Surf.val.j.i = NULL
    Surf.val = NULL}
  if(missing(Keep.other.grp)){Keep.other.grp = F}
  if(missing(Facet.size.egal)){Facet.size.egal = F}
  if(missing(Anomaly)){Anomaly = F}
  if(missing(Annotation)){Annotation = F}
  if(missing(Mono.core)){Mono.core = F}
  if(missing(Multi.clim)){Multi.clim = F}
  if(missing(Facet.scale.egal)){Facet.scale.egal = T}
  if(missing(Add.lim.space)){Add.lim.space = T}
  if(missing(Pollen.plot.merge)){Pollen.plot.merge = F}
  if(missing(Manual.y.val)){Manual.y.val = NULL}
  if(missing(Label.group)){Label.group = NULL}
  if(missing(Remove.models)){Remove.models = NULL}
  if(missing(Smooth.param)){Smooth.param = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Save.Rds)){Save.Rds = NULL}
  if(missing(Group.name)){Group.name = NULL}
  if(missing(Cores.lab)){Cores.lab = NULL}
  if(missing(Select.model)){Select.model = NULL}
  if(missing(Zone.clim)){Zone.clim = NULL}
  if(missing(Facet.T)){Facet.T = F}
  if(missing(Name.zone)){Name.zone = NULL}
  if(missing(Temp.zone)){Temp.zone = rep("C", length(Zone.clim))}
  if(missing(Clim.lab)){Clim.lab = NULL}
  if(missing(Age.select)){Age.select = "Age"}
  if(missing(Title.age)){Title.age = "Time (cal. year BP)"}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Select.interv)){Select.interv = 1000}
  if(missing(Zoom.box)){Zoom.box = NULL}
  if(missing(CRUTS)){CRUTS = NULL}
  if(missing(Pourc.clim.param)){Pourc.clim.param = NULL}
  if(missing(Extract.fitting)){Extract.fitting = F}
  if(missing(Smooth.SD)){Smooth.SD = T}
  
  #### Graphical settings ####
  Plot.list.col <- list()
  Cores <- names(MGDGT)
  Rect.color.scale <- c("W" = "#E76D51",
                        "C" = "#75AADB",
                        "G" = "grey40",
                        "D" = "#c67f05",
                        "Wt" = "#004266")
  
  Rect.color.scale <- Rect.color.scale[unique(Temp.zone)]
  
  if(Pollen.plot.merge == T){
    Clim.lab = NULL
    Title.just = NULL
  }
  else{Title.just = NULL}
  
  Legende.position.chart <- c("top", rep("none", (length(Cores)-1)))
  Legende.title.chart <- c(18, rep(0.5, (length(Cores)-1)))
  
  Legende.title.col.chart <- c("black", rep("white", (length(Cores)-1)))
  Legende.axis.chart <- c(rep("white", (length(Cores)-1)), "grey30")
  
  if(Mono.core == F){
    Legende.age.chart <- c(rep(0, (length(Cores)-1)), 11)
    Legende.Core.name.chart <- c(17, rep(0, (length(Pclim)-1)))
    Legende.biblio <- rep(3, length(Pclim))
    Axis.type.x <- element_blank()
    Axis.type.y <- element_blank()
  }
  else{
    Legende.position.chart <- c("top", rep("none", (length(Pclim)-1)))
    Axis.type.x <- element_line(colour = "grey30", lineend = "butt")
    Axis.type.y <- element_line(colour = "grey30", lineend = "butt")
    if(Multi.clim == F){Legende.age.chart <- c(rep(0, (length(Pclim)-1)), 11)}
    else{Legende.age.chart = 11}
    Legende.Core.name.chart = rep(18, length(Pclim))
    Legende.biblio <- c(2.5, rep(0, (length(Pclim)-1)))
  }
  
  if(is.null(Zone.clim) == F){
    yo = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                    xmax = Zone.clim[seq(2,length(Zone.clim), by=2)], 
                    Temp.col = Temp.zone)
    if(Zone.clim.dashed == T){LT <- 2}
    else{LT <- NA}
    Rect.clim <- geom_rect(data = yo, inherit.aes = F,
                           mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf, fill = Temp.col),
                           alpha = Temp.col.alpha, color = LT, linewidth = 0.3, linetype = 2)} 
  else{Rect.clim <- NULL}
  if(is.null(Label.group)==F){Label.cat <- Label.group[length(Label.group)]}
  else{Label.cat <- "A"}
  
  if(is.null(Name.zone) == F){yo2 = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                                               xmax = Zone.clim[seq(2,length(Zone.clim), by=2)],
                                               Temp.col = Temp.zone,
                                               Title.zone = Name.zone,
                                               Categorie = Label.cat
  )}
  else{yo2 = data.frame(xmin = 0, xmax = 0, Temp.col = "", Title.zone = "", Categorie = "")}
  
  if(is.null(Hiatus.age) == F){
    Hiatus.area <- data.frame(xmin = Hiatus.age[1], xmax = Hiatus.age[2])
    Hiatus.area <- geom_rect(data = Hiatus.area, inherit.aes = F,
                             mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill = "grey97",
                             alpha = 1, color = "grey60", linewidth = 0, linetype = 0)}
  else{Hiatus.area <- NULL}
  
  if(Show.x.axis == F){
    My_text_x <- element_blank()
    My_tick_x <- element_blank()
    Axis.type.x <- element_blank()
  }
  else{
    My_tick_x <- element_line(colour = "grey30")
    Axis.type.x <- element_line(colour = "grey30", lineend = "butt")
    Axis.type.y <- element_line(colour = "grey30", lineend = "butt")
  }
  
  if(Panel.box == T){
    My_panel.box <- element_rect(fill = NA, color = "grey30",)
    Axis.type.y <- element_blank()
  }
  else{My_panel.box <- element_blank()}
  
  #### Shape settings ####
  if(is.null(Shape.group) == T){myshapes <- c(16, 1, 3, 18, 13, 4, 15, 2, 19, 9, 7, 20)}
  else{myshapes <- Shape.group}
  names(myshapes) <- Label.group
  myshapes <- myshapes[1:length(Label.group)]
  if(is.null(Group.name) == F){names(Group.name) <- Label.group}
  
  #### Case only some Surface Values ####
  if(is.null(Surf.val) == F){A <- setNames(data.frame(matrix(ncol = length(setdiff(Pclim,names(Surf.val))), nrow = nrow(Surf.val))), setdiff(Pclim,names(Surf.val))) 
  row.names(A) <- row.names(Surf.val)
  A[is.na(A)] <- NA
  A <- cbind(Surf.val, A)
  Surf.val <- A[,sort(names(A))]}
  
  #### Convert shape of MGDGT matrix ####
  Convert.mat.gdgt <- function(M){
    M.reshape <- list()
    for(i in 1:length(M)){
      M.reshape[[i]] = sapply(Pclim, function(x) as.data.frame(M[[i]][grepl(x, names(M[[i]]))]), simplify = F)
      M.reshape[[i]] <- Map(cbind,Age = M[[i]][Age.select], M.reshape[[i]])
      names(M.reshape)[[i]] <- names(M)[i]
    }
    M.crop <- lapply(1:length(M.reshape[[1]]), function(i) lapply(M.reshape, "[[", i))
    names(M.crop) <- Pclim
    return(M.crop)
  }
  
  MGDGT.crop <- Convert.mat.gdgt(MGDGT)
  
  if(is.null(MMin) == F){MMin.crop <- Convert.mat.gdgt(MMin)}
  if(is.null(MMax) == F){MMax.crop <- Convert.mat.gdgt(MMax)}
  if(Extract.fitting == T){Select.models.tot <- list()}
  
  #### Main Loop ####
  if(is.null(Select.model) == F){
    for(a in 1:length(Select.model)){length(Select.model[[a]]) <- length(Cores)}
    Select.model <- data.frame(Select.model)}
  
  for(j in 1:length(Pclim)){
    Plot.list <- list()
    if(is.null(Surf.val) == F){Surf.val.j <- subset(Surf.val, select = c(Pclim[[j]]))}
    Mplot <- MGDGT.crop[[j]]
    
    if(Extract.fitting == T){Select.models <- data.frame(Age = seq(Limites[1], Limites[2], by = 10))}
    if(is.null(MMax) == F){Mplot.max <- MMax.crop[[j]]}
    if(is.null(MMin) == F){Mplot.min <- MMin.crop[[j]]}
    if(Mono.core == T){Legend.val <- Legende.age.chart[j]}
    for(i in 1:length(Cores)){
      #### Options ####
      Name.core <- Cores[i]
      Mplot.core <- Mplot[[i]]
      
      if(any(grepl(Pclim[j], names(Mplot.core))) == F){
        print(paste(Pclim[j], "not present in the GDGT matrix"))
        next}
      
      if(is.null(MMax) == F){Mplot.core.max <- Mplot.max[[i]]}
      if(is.null(MMin) == F){Mplot.core.min <- Mplot.min[[i]]}
      if(Mono.core == F){Legend.val <- Legende.age.chart[i]}
      if(is.null(Surf.val) == F){Surf.val.j.i <- Surf.val.j[Name.core,]}
      
      if(Show.x.axis == T){My_text_x <- element_text(angle = 45, hjust = 1, size = Legend.val, colour = "grey30")}
      
      #### CRUTS v4 plot ####
      if(is.null(CRUTS) == F){
        if(class(CRUTS) == "data.frame"){CRUTS.i <- CRUTS}
        else{
          CRUTS.i <- NULL
          CRUTS.i <- CRUTS[[Name.core]]
        }
        
        if(length(grep(Pclim[j], names(CRUTS.i))) > 0){
          CRUTS.i.j <- CRUTS.i[c(which(names(CRUTS.i) == Age.select), grep(Pclim[j], names(CRUTS.i)))]
          names(CRUTS.i.j)[2] <- "value"
          CRUTS.line <- geom_line(inherit.aes = F, data = CRUTS.i.j, mapping = aes(x = Age, y = value), color = "grey30", linetype = "solid", linewidth = 0.3, alpha = 0.6)
        }
        else{CRUTS.line <- NULL}
      }
      else{CRUTS.line <- NULL}
      
      #### Selection des modèles / Param.clim ####
      To.keep <- setdiff(names(Mplot.core), Remove.models)
      print(To.keep)
      if(is.null(Keep.models) == F){To.keep <- intersect(names(Mplot.core), c(Age.select, Keep.models))}
      Mplot.core <- Mplot.core[,To.keep]
      if(is.integer(Mplot.core)){next}
      
      if(is.null(MMin) == F){
        To.keep.min <- setdiff(names(Mplot.core.min), Remove.models)
        if(is.null(Keep.models) == F){To.keep.min <- intersect(names(Mplot.core.min), c(Age.select, Keep.models))}
        Mplot.core.min <- Mplot.core.min[,To.keep.min]
      }
      if(is.null(MMax) == F){Mplot.core.max <- Mplot.core.max[,To.keep.min]}
      if(is.null(Cores.lab) == F){Name.core <- Cores.lab}
      
      #### Specific selection of models ####
      Select.model.i <- Select.model[grep(Pclim[j], names(Select.model))]
      
      if(is.null(Select.model.i) == F){
        if(ncol(Select.model.i)>=1){
          To.keep <- Age.select
          for(k in 1:ncol(Select.model.i)){
            if(i %in% Select.model.i[,k] == T){
              To.keep <- c(To.keep, names(Select.model.i)[k])
            }
          }
        }
        Mplot.core <- Mplot.core[,To.keep]
      }
      
      #### Calcul en anomalies ####
      Lim.ano = NULL
      if(Anomaly == T & is.null(Surf.val) == F){
        Anomaly.conv.2 <- function(M){
          Keep.Age <- M[,1]
          M <- M[,-1]
          Msurf <- setNames(data.frame(matrix(ncol = ncol(M), nrow = nrow(M), Surf.val.j.i)), names(M))
          row.names(Msurf)<- row.names(M)
          M <- (Msurf - M)/Msurf
          M <- cbind(Age = Keep.Age, M)
          return(M)}
        
        Anomaly.conv <- function(M, Val.mean){
          if(missing(Val.mean)){Val.mean = NULL}
          
          Keep.Age <- M[,Age.select]
          Keep.names <- names(M)
          M <- data.frame(M[,!grepl(Age.select, names(M))])
          names(M) <- Keep.names[!grepl(Age.select, Keep.names)]
          
          if(is.null(Val.mean) == T){Ac.val <- unlist(M[c(1),])}
          else{
            Val.mean <- data.frame(Val.mean[,!grepl(Age.select, names(Val.mean))])
            
            Ac.val <- unlist(Val.mean[c(1),])
          }
          
          for(i in 1:ncol(M)){
            M[,i] <- M[,i] - Ac.val[i]
          }
          
          M <- cbind(Age = Keep.Age, M)
          return(M)}
        Save.not.anom <- Mplot.core
        if(is.null(MMin) == F){Mplot.core.min <- Anomaly.conv(Mplot.core.min, Val.mean = Save.not.anom)}
        if(is.null(MMax) == F){Mplot.core.max <- Anomaly.conv(Mplot.core.max, Val.mean = Save.not.anom)}
        Mplot.core <- Anomaly.conv(Mplot.core)
        if(is.null(CRUTS) == F){
          if(length(grep(Pclim[j], names(CRUTS))) > 0){
            CRUTS <- Anomaly.conv(CRUTS)}}
        
        Lim.ano <- c(min(Mplot.core[-1,]), max(Mplot.core[-1]))
        Surf.val.j.i <- 0
      }
      if(Anomaly == F & is.null(Surf.val) == T | Facet.scale.egal == F){Lim.ano = NULL}
      
      #### Correction surf line ####
      if(is.null(Surf.val) == F){
        if(is.na(Surf.val.j.i) == F){
          if(is.null(Surf.val.max.age) == T){
            Surf.line <- geom_hline(yintercept = Surf.val.j.i, linetype = "dotdash", linewidth = 0.55, color = "black", alpha = 0.55)
            Surf.dot <- NULL}
          else{
            if(Color.by.proxy == T){Col_surf <- "darkred"}
            else{Col_surf <- "#5E4FA2"}
            
            A = data.frame(x = Limites[1], xend = Surf.val.max.age, yend = Surf.val.j.i, y = Surf.val.j.i)
            Surf.line <- geom_segment(inherit.aes = F, data = A, aes(yend = yend, x = x, xend = xend, y = y), linetype = "dotdash", linewidth = 0.55, color = "black", alpha = 0.55)
            Surf.dot <- geom_point(inherit.aes = F, data = A, aes(x = x, y = y), color = Col_surf, alpha = 1, size = Dot.size*1.5, shape = 9)
          }
          
        }
        if(is.na(Surf.val.j.i) == T | is.null(Surf.val.j.i) == T){Surf.line <- NULL; Surf.dot <- NULL}}
      else{Surf.line <- NULL; Surf.dot <- NULL}
      
      #### Melt des données ####
      Mtot.melt <- melt(Mplot.core, id = "Age")
      
      if(is.null(MMin) == F){
        Mtot.melt.min <- melt(Mplot.core.min, id = "Age")
        Mtot.melt <- cbind(Mtot.melt, min.val = Mtot.melt.min$value)
      }
      else{Mtot.melt$min.val <- Mtot.melt$value}
      
      if(is.null(MMax) == F){
        Mtot.melt.max <- melt(Mplot.core.max, id = "Age")
        Mtot.melt <- cbind(Mtot.melt, max.val = Mtot.melt.max$value)
      }
      else{Mtot.melt$max.val <- Mtot.melt$value}
      
      #### Set the categories ####
      if(is.null(Label.group) == F){
        Categorie <- as.character(Mtot.melt$variable)
        for(k in Label.group){Categorie[grepl(k, Categorie)] <- k}
        
        if(Keep.other.grp == T){
          Add.to.other <- setdiff(unique(Categorie), Label.group)
          for(k in Add.to.other){Categorie[grepl(k, Categorie)] <- "Other"}
          Categorie <- factor(Categorie, levels = c(Label.group, "Other"))
          myshapes <- cbind(myshapes, Other = 20) # Dernier ajout, a vérif si ca marche avec les autres lake que Fazilman
        }
        else{Categorie <- factor(Categorie, levels = Label.group)}
        if(is.null(Group.name) == T){Group.name <- levels(Categorie)}
      }
      else{Categorie <- rep("A", nrow(Mtot.melt))}
      
      Mtot.melt <- cbind(Mtot.melt, Categorie)
      Mtot.melt <- na.omit(Mtot.melt)
      
      #### Facet param ####
      if(Facet.T == T){
        if(unique(yo2$Categorie) %in% unique(as.character(Mtot.melt$Categorie)) == F){
          yo2$Categorie <- unique(as.character(Mtot.melt$Categorie))}
        Facet.choix <- ~Categorie
      }
      else{
        Facet.choix <- NULL}
      
      if(Facet.scale.egal == T){Facet.scale <- NULL}
      else(Facet.scale <- "free_y")
      
      if(Facet.size.egal == T){
        if(length(setdiff(levels(Mtot.melt$Categorie), Mtot.melt$Categorie))>=1){
          Mtot.melt$variable <- factor(Mtot.melt$variable, levels = sort(levels(Mtot.melt$variable)))
          Mtot.melt$Categorie <- factor(Mtot.melt$Categorie, levels = unique(Mtot.melt$Categorie))
          Group.name <- levels(Mtot.melt$Categorie)
        }
      }
      #### Limites Settings ####
      if(missing(Limites)){Limites = c(round(min(Mtot.melt$Age),digits = -2), max(Mtot.melt$Age))}
      if(Add.lim.space == F){Limites.zone <- c(Limites[1], Limites[2])}
      else{Limites.zone <- c(Limites[1], (Limites[2]+0.23*Limites[2])) }
      
      Mtot.melt <- Mtot.melt[Mtot.melt$Age <= Limites[2] & Mtot.melt$Age >= Limites[1],]
      
      #### Set the etiquettes des lignes ####
      Lab.etiquette <- rep("", nrow(Mtot.melt))
      Mtot.melt <- cbind(Mtot.melt, Lab.etiquette)
      Max.age <- max(Mtot.melt$Age)
      Ind.max <- which(Mtot.melt$Age %in% Max.age)
      Mtot.melt[Ind.max, ncol(Mtot.melt)] <- as.character(Mtot.melt[Ind.max, 2])
      
      if(length(grep("_", Mtot.melt$Lab.etiquette))>0){
        Mtot.melt$Lab.etiquette = paste(sub("_", "[", Mtot.melt$Lab.etiquette), "]", sep = "")
        Mtot.melt$Lab.etiquette = sub("_", "~", Mtot.melt$Lab.etiquette)
        Mtot.melt$Lab.etiquette = sub("_", "-", Mtot.melt$Lab.etiquette)
      }
      
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "MAAT[DJ~5Me]"] <- "MAAT[DJ~5*Me]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "CBT5Me"] <- "CBT[5][Me]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "CBT6Me"] <- "CBT[6][Me]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "CBTp"] <- "CBT*minute"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "MBTp"] <- "MBT*minute"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "CBT5MeDJ"] <- "CBT[5][MeDJ]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "MBTp5MeDJ"] <- "MBT*minute[5][Me~DJ]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "MBTp5Me"] <- "MBT*minute[5][Me]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "MBTp6Me"] <- "MBT*minute[6][Me]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "CBT.IR"] <- "IR[6][Me]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "CBT.Ib.frac"] <- "f(Ib)"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "CBT.Ic.frac"] <- "f(Ic)"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "CBT.IIa.frac"] <- "f(IIa)"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "Rib"] <- "R[i/b]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "Rib.CI"] <- "CI"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "Rib.CI.community"] <- "CI"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "Rib.Index2"] <- "Index[2]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "Rib.Index2.community"] <- "Index[2]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "Rib.GDGT0.Crenar"] <- "GDGT[0]/Cren"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "Tetra.methyl"] <- "\'%\'[tetra]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "Penta.methyl"] <- "\'%\'[penta]"
      Mtot.melt$Lab.etiquette[Mtot.melt$Lab.etiquette == "Hexa.methyl"] <- "\'%\'[hexa]"
      
      Mtot.melt$Lab.etiquette[str_count(Mtot.melt$Lab.etiquette) <= 1] <- NA
      
      #### Color settings ####
      if(Color.by.proxy == F){
        A = 1:11
        Keep.col <- A[-c(6,7)] #there are 9, I exluded the two lighter hues
        if(length(levels(Mtot.melt$variable)) <= 3){Keep.col <- A[-c(5,6,7,8)]}
        if(length(levels(Mtot.melt$variable)) == 1){Keep.col <- A[-c(1,2,3,4,5,6,7,8,9,10)]}
        if(length(levels(Mtot.melt$variable)) == 2){Keep.col <- A[-c(1,3,4,5,6,7,8,9,10)]}
        
        my_orange = brewer.pal(n = 11, "Spectral")[Keep.col] 
        orange_palette = colorRampPalette(my_orange)
        my_orange = orange_palette(length(unique(Mtot.melt$variable)))
        
        if(is.null(Manual.color.scale) == F){my_orange <- Manual.color.scale}
        
        My_courbe <- geom_line(data = Mtot.melt, mapping = aes(color = variable), linetype = "solid", linewidth = 0.5, alpha = 0.6)
        My_dots <- geom_point(aes(color = variable), size = Dot.size)
      }
      else{
        my_orange <- "#CF5612"
        My_courbe <- geom_line(data = Mtot.melt, color = my_orange, linetype = "solid", linewidth = 0.5, alpha = 0.6)
        My_dots <- geom_point(color = my_orange, size = Dot.size)
      }
      
      
      #### Mono.core settings ####
      if(Mono.core == T & is.null(Clim.lab) == F){
        Legende.core.clim <- Legende.position.chart[j]
        Name.core <- Clim.lab[j] # Bug possible ici pour les noms.
        Clim.lab[[j]] <- ""
        if(j>1){yo2$Title.zone <- ""}
      }
      if(Mono.core == T & is.null(Clim.lab) == T){Legende.core.clim <- Legende.position.chart[j]}
      if(Mono.core == T){
        Break.special <- scale_y_continuous(name = Name.core[[i]],  breaks = scales::pretty_breaks(n = 4), limits = NULL)
        Axis.special <- NULL}
      if(Mono.core == F){
        if(is.null(Manual.y.val) == T){
          #### Breaks definition ####
          Kmax <- max(Mtot.melt$value, Surf.val.j.i, na.rm = T)
          Kmin <- min(Mtot.melt$value, Surf.val.j.i, na.rm = T)
          if(Kmax >= 40){
            Kmin <- round(Kmin, digits = -1)
            Kmax <- round(Kmax, digits = -1)}
          else{
            Kmin <- round(Kmin, digits = 1)
            Kmax <- round(Kmax, digits = 1)
          }
          Kmid <- Kmin + (Kmax - Kmin)*1/3
          Kmid2 <- Kmin + (Kmax - Kmin)*2/3
          Range.round <- c(Kmin, Kmax)
          
          if(Kmax >= 40){
            Range.round <- round(Range.round, digits = -1)
          }
          else{
            if(abs(Kmax-Kmin)<10){Range.round <- round(Range.round, digits = 0)}
            else{Range.round <- round(Range.round, digits = -1)}
          }}
        else{
          Clim.break <- gsub("\\..*", "", names(Manual.y.val))
          Core.break <- gsub(".*\\.", "", names(Manual.y.val))
          Match.clim <- which(Clim.break %in% Pclim[j])
          Match.core <- which(Core.break %in% Cores[i])
          Match.tot <- unique(intersect(Match.clim, Match.core), intersect(Match.core, Match.clim))
          Range.round <- Manual.y.val[[Match.tot]]
        }
        
        #### Axes ####
        if(Show.x.axis == F){My_data <- data.frame(A = Range.round[c(1, length(Range.round))], B = c(0,0))}
        else{My_data <- data.frame(A = Range.round[c(1, length(Range.round))], B = Limites)}
        
        Axis.special <- geom_rangeframe(data = My_data, 
                                        inherit.aes = F, mapping = aes(x = B, y = A), 
                                        colour = "grey50", sides = "lb", size = 0.5) 
        Legende.core.clim <- Legende.position.chart[i]
        
        if(length(Range.round) > 2){
          Break.special <- scale_y_continuous(name = Name.core[[i]], breaks = Range.round, limits = NULL)
        }
        else{
          Break.special <- scale_y_continuous(name = Name.core[[i]], breaks = round(seq(Range.round[1], Range.round[2], length.out = 5), digits = 0) , limits = NULL)}
      }
      if(Pollen.plot.merge == T & Mono.core == T){Name.core <- paste(Name.core, " (brGDGT)", sep = "")}
      if(is.null(Smooth.param) == F){
        if(Color.by.proxy == F){
          Smooth.line <- geom_smooth(method = "loess", se = Smooth.SD, formula = "y ~ x", fullrange = F, level = 0.95, linetype="solid",
                                     linewidth = .8, show.legend = F, aes(fill = variable, colour = variable), span = Smooth.param, alpha = 0.2)}
        else{Smooth.line <- geom_smooth(method = "loess", se = Smooth.SD, formula = "y ~ x", fullrange = F, level = 0.95, linetype="solid",
                                        linewidth = .8, show.legend = F, fill = my_orange, color = my_orange, span = Smooth.param, alpha = 0.2)}
        
      }
      else{Smooth.line <- NULL}
      
      if(length(na.omit(unique(Mtot.melt$Lab.etiquette))) == 1){NY = 10}
      else{NY = 0}
      
      if(is.null(Display.legends) == F){Legende.core.clim <- Display.legends}
      #### Repel settings ####
      if(Show.repels == T){
        if(is.null(Repel.x) == T){Repel.x <- c(Limites[2], Limites.zone[2])}
        
        My_etiquettes_repel <-  geom_text_repel(mapping = aes(x = Age, label = Lab.etiquette, colour = variable), #nudge_y = NY, # ETIQUETTE MODELS
                                                direction = "y", hjust = 1, segment.square = F, segment.shape = 1, segment.linetype = 4, 
                                                force = Repel.repoussage, na.rm = T, xlim = Repel.x,
                                                segment.curvature = 0.2, segment.inflect = T, min.segment.length = 0.5,
                                                size = 4.5, parse = T, segment.size = 0.1, segment.colour = "grey50"
        )
      }
      else{My_etiquettes_repel <- NULL}
      
      #### RMSE barres ####
      if(is.null(RMSE.barres) == F){
        Last <- Mtot.melt[!is.na(Mtot.melt$Lab.etiquette), c("Age", "variable", "value", "Categorie")]
        Last$RMSE <- RMSE.barres$RMSE[match(Last$variable, row.names(RMSE.barres))]
        Last$max.val <- Last$value + Last$RMSE
        Last$min.val <- Last$value - Last$RMSE
        
        Last$Age <- Limites.zone[2]
        Last$epsilon <- 0.01*Last$Age
        if(nrow(Last) > 1){for(ii in 2:nrow(Last)){Last$Age[ii] <- Last$Age[ii-1] - Last$epsilon[ii-1]}}
        
        Seg.RMSE <- geom_segment(inherit.aes = F, data = Last, mapping = aes(x = Age, xend = Age, y = max.val, yend = min.val, color = variable), linetype = "solid", linewidth = 1, alpha = 1)
      }
      else{Seg.RMSE <- NULL}
      #### Zoom area ####
      if(is.null(Zoom.box) == F){
        
        Select.in.zone <- Mtot.melt[Mtot.melt$Age >= Zoom.box[1] & Mtot.melt$Age <= Zoom.box[2],]
        
        Zoom.area <- data.frame(xmin = Zoom.box[1], xmax = Zoom.box[2], ymax = max(Select.in.zone$max.val, na.rm = T), ymin = min(Select.in.zone$min.val, na.rm = T))
        Zoom.area <- geom_rect(data = Zoom.area, inherit.aes = F,
                               mapping = aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax), fill = NA,
                               alpha = 1, color = "grey30", linewidth = .8, linetype = 2)}
      else{Zoom.area <- NULL}
      
      #### PLOTS ####
      Plot.list[[i]] <- ggplot(data = Mtot.melt, 
                               mapping = aes(x = Age, y = value, shape = as.factor(Categorie)), fill = variable)+#, color = DB)) +
        #### Rectangles ####
      scale_fill_manual(values = Rect.color.scale, guide = "none")+
        Rect.clim +
        geom_vline(xintercept = Manual.vlines, col = "grey30", lty = 2, alpha = 0.7)+
        
        #### Items ####
      new_scale_fill()+
        scale_shape_manual(values = myshapes, label = Group.name) +
        scale_color_manual(values = my_orange, guide = "none")+
        scale_fill_manual(values = my_orange, guide = "none")+
        Smooth.line + CRUTS.line +
        Seg.RMSE +
        My_courbe +     # dashted, dotted, solid
        geom_ribbon(mapping = aes(ymin = min.val, ymax = max.val, fill = variable), alpha = 0.1, size = 0.15, linetype = "dashed") +
        new_scale_fill()+ Hiatus.area +
        My_dots +
        scale_x_continuous(name = Title.age, breaks = c(Limites[1], round(seq(0, Limites[2], by = Select.interv))), expand = c(0.02,0.02))+
        coord_cartesian(xlim = Limites.zone, ylim = Lim.ano, clip = "off")+
        guides(shape = guide_legend(override.aes = list(size = Legende.biblio[j], alpha = 0.5)))+
        facet_wrap(Facet.choix, ncol = 1, scales = Facet.scale, drop = F)+
        ggtitle(Clim.lab[[j]])+
        Zoom.area +
        
        #### Annotations ####
      My_etiquettes_repel +
        Surf.line + Surf.dot +
        new_scale_color()+
        scale_color_manual(values = Rect.color.scale, guide = "none", name = NULL, labels = NULL, breaks = NULL, na.translate = FALSE)+
        geom_text(data = yo2, inherit.aes = F,   # PERIODE CLIMATIQUE TITRE 
                  aes(x = (xmax+xmin)/2, y = Inf, label = Title.zone, color = Temp.col), size = 3.5, vjust = 1.5, angle = Lab.clim.angle, fontface = "bold") +
        Axis.special +
        #### Theme ####
      Break.special +
        theme(
          plot.title =  element_text(size = Legende.title.chart[i], hjust = 0.5, color = Legende.title.col.chart[i]),
          plot.background = element_blank(),
          axis.text.x = My_text_x,
          axis.text.y = element_text(size = 11, colour = "grey30"),
          axis.title.x = element_text(size = Legend.val, colour = "grey20"),           # taile police du titre de l'axe
          axis.title.y = element_text(size = Legende.Core.name.chart[j], vjust = Title.just, colour = "grey20"),     # taile police du titre de l'axe
          axis.line.x = Axis.type.x, 
          axis.line.y = Axis.type.y,
          
          axis.ticks.x = My_tick_x,
          axis.ticks.y = element_line(colour = "grey30"),
          legend.title = element_blank(),
          legend.key = element_blank(),  # carré autours du symbole
          legend.position = Legende.core.clim,            # Legendes DB
          legend.justification = c("left"),               # left, top, right, bottom
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.text = element_text(size = Legende.biblio[j]*4, color = "grey20", hjust = 0),
          panel.background = My_panel.box,
          panel.spacing.x = unit(-3, "cm"),
          panel.spacing.y = unit(0, "cm"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          strip.text.x = element_blank(),
          strip.placement = "none",
          strip.background = element_blank()
        )
      
      #### Extract fitting ####
      if(Extract.fitting == T){
        if(length(unique(Mtot.melt$variable)) == 1){
          #### New version (loess) ####
          age_grid <- seq(Limites[1], Limites[2], by = 10)
          loess_fit <- loess(formula(paste("value~", Age.select)), data = Mtot.melt, span = Smooth.param)
          Save.fit <- data.frame(age_grid, predict(loess_fit, newdata = data.frame(Age = age_grid)))
          if(Detrending.fitting == T){Save.fit[[2]] <- resid(lm(Save.fit[[2]] ~ Save.fit[[1]], na.action = na.exclude))}
          Save.fit[[2]] <- 2*(Save.fit[[2]]-min(Save.fit[[2]], na.rm = T))/(max(Save.fit[[2]], na.rm = T)-min(Save.fit[[2]], na.rm = T))-1
          colnames(Save.fit) <- c(Age.select, paste(Cores[[i]], Pclim[[j]], unique(Mtot.melt$Model), sep = "."))
          
          #### Export ####
          Merge.mat <- merge(Select.models, Save.fit, by ="Age", all.x = T)
          Select.models <- Merge.mat
        }
        else(print("To calculate the climatic phase, please only select one model for each cores / param. clim."))
      }
    }
    #### Facet plot clim / cores ####
    if(Mono.core == T){Col.row <- " + "}
    else{Col.row <- " / "}
    Formula.patch <- paste(paste("Plot.list[[", seq(length(Cores)), "]]", sep = ""), collapse = Col.row)
    Plot.list.col[[j]] <- eval(parse(text = Formula.patch))
    if(Extract.fitting == T){Select.models.tot[[j]] <- Select.models}
  }
  
  #### Save plots and export ####
  if(Mono.core == T & Multi.clim == F){Col.row <- " / "}
  else{Col.row <- " | "}
  Formula.patch <- paste(paste("Plot.list.col[[", seq(length(Pclim)), "]]", sep = ""), collapse = Col.row)
  Ptot <- eval(parse(text = Formula.patch))
  
  if(Annotation == T){Ptot <- Ptot + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 18, vjust = 1.3, hjust = -0.3))}
  
  if(Mono.core == T & Multi.clim == F){
    Ptot <- Ptot + plot_layout(heights = Pourc.clim.param)
    Marges.normal <- c(-0.5,0.3,0,0) # top, right, bottom, left
    Ptot <- Ptot & theme(plot.margin = unit(Marges.normal, "cm"))}
  else{Ptot <- Ptot + plot_layout(widths = Pourc.clim.param)}
  
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(Ptot, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.plot)}}
  
  if(is.null(Save.Rds) == F & Extract.fitting == T){saveRDS(Select.models.tot, Save.Rds)}
  
  return(Ptot)
  
}

Plot.FT.summary <- function(FT, Cores, Pclim, Model, Select.interv, Cores.lab, Facet.T, Repel.T, Zoom.box = NULL, Legend.pos = NULL, RMSE.barres = NULL,
                            Temp.col.alpha = 0.1, Surf.val, Anomaly, GDGT.plot.merge, Smooth.T = T, Axis.right = F, Label.group, Smooth.param, Phase.clim, Title.display = T,
                            Surf.val.max.age = NULL, Detrending.fitting = F,
                            Save.Rds, Save.plot, W, H, Limites, Add.lim.space, Mono.core, Condensed, Mean.models, Save.plot.Rds = NULL, Dot.size = 2, Name.zone.angle = 0,
                            Manual.vlines = NULL, Smooth.sd, Only.fit, Zone.clim, Zone.clim.box = T, Name.zone, Merge.legends = F, Temp.zone, Clim.lab, Manual.y.val){
  #### Init Val ####
  if(missing(FT)){warning("Import the transfer function for plotting.")}
  if(missing(Cores)){warning("Select and sort the core(s) to plot.")}
  if(missing(Pclim)){warning("Select and sort the climate parameter(s) to plot.")}
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(ggthemes)
  library(ggrepel) # nom des lignes à côté
  
  if(missing(Surf.val)){
    Surf.val.j.i = NULL
    Surf.val = NULL}
  if(missing(Mean.models)){Mean.models = F}
  if(missing(Anomaly)){Anomaly = F}
  if(missing(Condensed)){Condensed = F}
  if(missing(Mono.core)){Mono.core = F}
  if(missing(Add.lim.space)){Add.lim.space = T}
  if(missing(Manual.y.val)){Manual.y.val = NULL}
  if(missing(Cores.lab)){Cores.lab = NULL}
  if(missing(GDGT.plot.merge)){GDGT.plot.merge = F}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Save.Rds)){Save.Rds = NULL}
  if(missing(Label.group)){
    Label.group <- c("NMSDB", "MDB", "COST", "COSTDB", "EAPDB", "TUDB", "WASTDB", "STDB", "ACADB")}
  if(missing(Zone.clim)){Zone.clim = NULL}
  if(missing(Name.zone)){Name.zone = NULL}
  if(missing(Temp.zone)){Temp.zone = rep("C", length(Zone.clim))}
  if(missing(Clim.lab)){Clim.lab = NULL}
  if(missing(Limites)){Limites = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Select.interv)){Select.interv = 1000}
  if(missing(Smooth.param)){Smooth.param = 0.6}
  if(missing(Only.fit)){Only.fit = F}
  if(missing(Facet.T)){Facet.T = T}
  if(missing(Repel.T)){Repel.T = F}
  if(missing(Phase.clim)){Phase.clim = F}
  if(missing(Smooth.sd)){Smooth.sd = T}
  if(missing(Model)){Model = c("BRT", "MAT", "WAPLS")}
  
  #### Graphical settings ####
  Plot.list.col <- list()
  Rect.color.scale <- c("W" = "#E76D51",
                        "C" = "#75AADB",
                        "G" = "grey40",
                        "D" = "#c67f05",
                        "Wt" = "#004266")
  Rect.color.scale <- Rect.color.scale[unique(Temp.zone)]
  
  if(GDGT.plot.merge == T){
    Title.age = NULL
    Legende.age.chart <- rep(0, length(Cores))
    if(Add.lim.space == T){
      FT <- lapply(FT, function(x) x[which(x$Age <= Limites[2]),])
      Limites[2] <- (Limites[2] + Limites[2]*0.23)}
    Legende.Core.name.chart <- c(17, rep(0, (length(Pclim)-1)))
  }
  if(GDGT.plot.merge == F){
    Legende.age.chart <- c(rep(0, (length(Cores)-1)), 11)
    Title.age = "Time (cal. year BP)"
    Legende.Core.name.chart <- c(16, rep(0, (length(Pclim)-1)))
  }
  Legende.position.chart <- c("top", rep("none", (length(Cores)-1)))
  if(is.null(Legend.pos) == F){Legende.position.chart <- Legend.pos}
  
  Legende.title.chart <- c(16, rep(0.5, (length(Cores)-1)))
  
  Legende.title.col.chart <- c("black", rep("white", (length(Cores)-1)))
  Legende.axis.chart <- c(rep("white", (length(Cores)-1)), "grey")
  if(length(Smooth.param) != length(Cores)){Smooth.param <- rep(Smooth.param, length(Cores)/length(Smooth.param))}
  if(is.null(Zone.clim) ==F){
    if(Zone.clim.box == T){Box.col <- "grey"}
    else{Box.col <- NA}
    yo = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                    xmax = Zone.clim[seq(2,length(Zone.clim), by=2)], 
                    Temp.col = Temp.zone)
    Rect.clim <- geom_rect(data = yo, inherit.aes = F,
                           mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf, fill = Temp.col),
                           alpha = Temp.col.alpha, color = Box.col, size = 0.3, linetype = 2, na.rm = T)} 
  else{Rect.clim <- NULL}
  if(is.null(Name.zone) == F){yo2 = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                                               xmax = Zone.clim[seq(2,length(Zone.clim), by=2)],
                                               Temp.col = Temp.zone,
                                               Title.zone = Name.zone,
                                               Model = Model[length(Model)])}
  else{yo2 = data.frame(xmin = 0, xmax = 0, Temp.col = "", Title.zone = "")}
  
  
  
  if(Facet.T == T){
    Geom.point <- geom_point(shape = 20, size = Dot.size*.9)
    Facet.disp <- facet_grid(rows = vars(Model), switch = "y")}
  else{
    Geom.point <- geom_point(size = Dot.size)
    Facet.disp <- facet_grid(rows = NULL, switch = "y")}
  
  
  #### Case only some Surface Values ####
  if(is.null(Surf.val) == F){A <- setNames(data.frame(matrix(ncol = length(setdiff(Pclim,names(Surf.val))), nrow = nrow(Surf.val))), setdiff(Pclim,names(Surf.val))) 
  row.names(A) <- row.names(Surf.val)
  A[is.na(A)] <- NA
  A <- cbind(Surf.val, A)
  Surf.val <- A[,sort(names(A))]}
  
  #### Reorganized the surface values according to cores ####
  
  #### Main Loop ####
  FT.crop <-sapply(Cores, function(x) as.data.frame(FT[grepl(x, names(FT))]), simplify = F)
  
  if(is.null(RMSE.barres) == F){
    RMSE.barres <- RMSE.barres[grepl("Param", names(RMSE.barres))]
    KRN <- row.names(RMSE.barres[[1]])
    RMSE.barres <- data.frame(sapply(RMSE.barres, function(x) x$RMSE.choice))
    row.names(RMSE.barres) <- KRN
    names(RMSE.barres) <- gsub("\\.Best\\.Param", "", names(RMSE.barres))
  }
  
  if(Phase.clim == T){Select.models.tot <- list()}
  for(j in 1:length(Pclim)){
    #### Settle surf val j ####
    print(paste("Ploting the", Pclim[j]))
    Plot.list <- list()
    if(is.null(Surf.val) == F){Surf.val.j <- subset(Surf.val, select = c(Pclim[j]))}
    if(Phase.clim == T){
      if(is.null(Limites) == T){stop("The limites have to be settled to apply Phase.clim = T")}
      if(Limites[2] > 100){Select.models <- data.frame(Age = seq(Limites[1], Limites[2], by = 10))}
      if(Limites[2] <= 100){Select.models <- data.frame(Age = seq(Limites[1], Limites[2], by = .01))}
    }
    
    if(is.null(RMSE.barres) == F){RMSE.barres.j <- RMSE.barres[row.names(RMSE.barres) == Pclim[j],]}
    
    #### Verbose ####
    library(lubridate)
    if(length(Cores) > 1){
      pb = txtProgressBar(min = 1, 
                          max = length(Cores), 
                          width = 40,
                          initial = 0,  style = 3) 
      
      init <- numeric(length(Cores))
      end <- numeric(length(Cores))}
    
    #### Loop 2 ####
    for(i in 1:length(Cores)){
      #### Settle surf val i ####
      if(length(Cores) > 1){init[i] <- Sys.time()}
      Name.core <- Cores[i]
      if(is.null(Surf.val) == F){
        Surf.val.j.i <- Surf.val.j[Name.core,]
        if(is.na(Surf.val.j.i) == T){Surf.val.j.i = NULL}
      }
      
      if(is.null(Cores.lab) == F){Name.core <- Cores.lab}
      
      #### Fusion des list en 1 unique matrice par Carotte ####
      Mtest <- FT.crop[[i]]
      Mtest <- Mtest[,!grepl("SEP", colnames(Mtest))]
      colnames(Mtest)[1]<-"Age"
      Mtest <- Mtest[,!grepl("\\.Age", colnames(Mtest))]
      
      #### Selection des modèles / Param.clim ####
      names(Mtest) <- sub("Psum", "SUMMERPR", names(Mtest)) 
      Ttete <- unlist(c("Age", sapply(Model, function(x) names(Mtest)[grepl(x, names(Mtest))])))
      Ttete <- c("Age", sapply(Pclim[[j]], function(x) Ttete[grepl(x, Ttete)]))
      Mtest <- Mtest[Ttete]
      
      #### Calcul en anomalies ####
      if(Anomaly == T & is.null(Surf.val) == F){
        Keep.Age <- Mtest[,1]
        Msurf <- setNames(data.frame(matrix(ncol = ncol(Mtest), nrow = nrow(Mtest), Surf.val.j.i)), names(Mtest))
        row.names(Msurf)<- row.names(Mtest)
        Mtest <- (Mtest - Msurf)/Msurf
        Surf.val.j.i <- 0
        Lim.ano <- c(min(Mtest[2], na.rm = T), max(Mtest[2], na.rm = T))
        Lim.ano <- c(min(Lim.ano[1], Surf.val.j.i, na.rm = T), max(Lim.ano[2], Surf.val.j.i, na.rm = T))
        Mtest <- cbind(Age = Keep.Age, Mtest)
      }
      else{Lim.ano = NULL}
      
      #### Melt des données ####
      Mtot.melt <- melt(Mtest, id = "Age")
      Categorie <- t(data.frame(strsplit(as.character(Mtot.melt$variable),split="\\.")))
      Mtot.melt <- cbind(Mtot.melt, Categorie)
      Mtot.melt <- Mtot.melt[,-2]
      colnames(Mtot.melt) <- c("Age", "variable", "Lake", "Model", "DB", "Param.clim")
      Mtot.melt$DB <- factor(Mtot.melt$DB, levels = Label.group)
      
      row.names(Mtot.melt) <- paste("P", 1:nrow(Mtot.melt), sep = "")
      Mtot.melt <- na.omit(Mtot.melt)
      
      #### Last graphical settings ####
      if(Condensed == T){
        if(i < length(Cores)){
          Xtick <- element_blank()
          Xline <- element_blank()}
        else{
          Xtick <- element_line(colour = "grey55")
          Xline <- element_line(colour = "grey55", lineend = "butt")}
        if(i > 1){yo2$Title.zone.i <- ""}
        else{yo2$Title.zone.i <- yo2$Title.zone}
      }
      else{
        Xtick <- element_line(colour = "grey55")
        Xline <- element_line(colour = "grey55", lineend = "butt")
        yo2$Title.zone.i <- yo2$Title.zone}
      if(GDGT.plot.merge == T){
        Xtick <- element_blank()
        Xline <- element_blank()
        if(is.null(Cores.lab) == T){
          if(Mono.core == T){Name.core <- "Pollen-based \n climate reconstructions"}
          else{Name.core <- paste(Name.core, " (pollen)", sep = "")}}
      }
      
      if(length(Name.core) == 1){Name.core.plot <- Name.core}
      else{Name.core.plot <- Name.core[[i]]}
      
      if(Title.display == F){My_title <- element_blank()}
      else{My_title <- element_text(size = Legende.title.chart[i], hjust = 0.5, color = Legende.title.col.chart[i])}
      #### Color settings ####
      My_color <- c("WASTDB" = "#e2a064ff", "WAST" = "#e2a064ff", 
                    "NMSDB" = "#F3A481", "MDB" = "#963326", "TUDB" = "#0094AF", "TUSDB" = "#0094AF",
                    "ST" = "#91C4DD", "STDB" = "#91C4DD", 
                    "ACADB" = "#963326", "CAUCDB" = "#0F3361",
                    "COSTDB" = "#f3c768ff", "COST" = "#f3c768ff",
                    "MEDTEMP" = "#F3A481", "TEMPSCAND" = "#91C4DD", 
                    "EAPDB" = "#0F3361", "TAIGDB" = "#32156eff")
      #### Limites ####
      if(is.null(Limites) == T){Limites = c(min(Mtot.melt$Age, na.rm = T), max(Mtot.melt$Age, na.rm = T))}
      else{
        Mtot.melt <- Mtot.melt[Mtot.melt$Age <= Limites[2],]
        Mtot.melt <- Mtot.melt[Mtot.melt$Age >= Limites[1],]
      }
      
      
      #### Add RMSE barres ####
      if(is.null(RMSE.barres) == F){
        RMSE.barres.j <- data.frame(t(RMSE.barres.j))
        # RMSE.barres.j$Age <- Limites[2]
        RMSE.barres.j$Age <- max(Mtot.melt$Age, na.rm = T)
        
        RMSE.barres.j$Lake <- Cores[i]
        RMSE.barres.j$Model <- gsub(".*\\.", "", row.names(RMSE.barres.j))
        RMSE.barres.j$DB <- gsub("\\..*", "", row.names(RMSE.barres.j))
        RMSE.barres.j$Param.clim <- names(RMSE.barres.j)[1]
        names(RMSE.barres.j)[1] <- "variable"
        rows_A <- apply(RMSE.barres.j[c(2:6)], 1, paste, collapse = "_")
        rows_B <- apply(Mtot.melt[c(1,3:6)], 1, paste, collapse = "_")
        matching_indices <- match(rows_A, rows_B)
        
        RMSE.barres.j$value <- Mtot.melt$variable[matching_indices]
        RMSE.barres.j$max.val <- RMSE.barres.j$value + RMSE.barres.j$variable
        RMSE.barres.j$min.val <- RMSE.barres.j$value - RMSE.barres.j$variable
        
        RMSE.barres.j <- RMSE.barres.j[!is.na(RMSE.barres.j$value),]
        
        RMSE.barres.j$epsilon <- 0.015*RMSE.barres.j$Age
        
        for(k in unique(RMSE.barres.j$Model)){
          TT <- RMSE.barres.j[RMSE.barres.j$Model == k,]
          for(kk in 1:nrow(TT)){
            TT[kk,]$Age <- TT[kk,]$Age + kk* TT[kk,]$epsilon   
          }
          RMSE.barres.j[RMSE.barres.j$Model == k,] <- TT
        }
        Seg.RMSE <- geom_segment(inherit.aes = F, data = RMSE.barres.j, mapping = aes(x = Age, xend = Age, y = max.val, yend = min.val, color = DB), linetype = "solid", linewidth = 1, alpha = 1)
        
      }
      else{Seg.RMSE <- NULL}
      
      #### Repels ####
      if(Repel.T == T){
        Mtot.melt$Lab <- NA
        Mtot.melt$Lab[which(Mtot.melt$Age %in% max(Mtot.melt$Age, na.rm = T))] <- paste(as.character(Mtot.melt$Param.clim[which(Mtot.melt$Age %in% max(Mtot.melt$Age, na.rm = T))]), "[",
                                                                                        as.character(Mtot.melt$Model[which(Mtot.melt$Age %in% max(Mtot.melt$Age, na.rm = T))]), "-",
                                                                                        as.character(Mtot.melt$DB[which(Mtot.melt$Age %in% max(Mtot.melt$Age, na.rm = T))]), "]", sep ="")
        if(length(na.omit(unique(Mtot.melt$Lab))) == 1){NY = 10}
        else{NY = 0}
        Repel <-  geom_text_repel(mapping = aes(x = Age, label = Lab),  nudge_y = NY, # ETIQUETTE MODELS
                                  force = 4, nudge_x  = 20000, direction = "y", hjust = 1,
                                  size = 4.5, parse = T, segment.size = 0.18, segment.colour = "grey70")
      }
      else{Repel <- NULL}
      
      #### Breaks ####
      if(is.null(Manual.y.val) == T){
        Kmax <- max(Mtot.melt$variable, Surf.val.j.i, na.rm = T)
        Kmin <- min(Mtot.melt$variable, Surf.val.j.i, na.rm = T)
        if(Kmax >= 40){
          Kmin <- round(Kmin, digits = -1)
          Kmin = floor(Kmin / 50) * 50
          Kmax <- round(Kmax, digits = -1)
          Kmax = ceiling(Kmax / 50) * 50
          Kmid <- Kmin + (Kmax - Kmin)*1/3
          Kmid = ceiling(Kmid / 50) * 50
          Kmid2 <- Kmin + (Kmax - Kmin)*2/3
          Kmid2 = ceiling(Kmid2 / 50) * 50
        }
        if(Anomaly == T & is.null(Surf.val) == F){
          Kmid <- Kmin + (Kmax - Kmin)*1/3
          Kmid2 <- Kmin + (Kmax - Kmin)*2/3
        }
        else{
          Kmin <- round(Kmin, digits = 1)
          Kmin = floor(Kmin / 0.5) * 0.5
          Kmax <- round(Kmax, digits = 1)
          Kmax = ceiling(Kmax / 0.5) * 0.5
          Kmid <- Kmin + (Kmax - Kmin)*1/3
          Kmid = ceiling(Kmid / 0.5) * 0.5
          Kmid2 <- Kmin + (Kmax - Kmin)*2/3
          Kmid2 = ceiling(Kmid2 / 0.5) * 0.5
        }
        Range.round <- c(Kmin, Kmid, Kmid2, Kmax)
        if(Kmax >= 40){
          Range.round <- round(Range.round, digits = -1)
        }
        if(Anomaly == T & is.null(Surf.val) == F){
          Range.round <- round(Range.round, digits = 2)
        }
        else{
          Range.round <- round(Range.round, digits = 1)
        }
        
      }
      else{
        Clim.break <- gsub("\\..*", "", names(Manual.y.val))
        Core.break <- gsub(".*\\.", "", names(Manual.y.val))
        Match.clim <- which(Clim.break %in% Pclim[j])
        Match.core <- which(Core.break %in% Cores[i])
        Match.tot <- unique(intersect(Match.clim, Match.core), intersect(Match.core, Match.clim))
        Range.round <- Manual.y.val[[Match.tot]]
      }
      
      My_axis_clim <- data.frame(A = Range.round[c(1,length(Range.round))], B = Limites)
      if(Axis.right == F){Axis.side <- "l"; Axis.side2 <- "left"}
      else{Axis.side <- "r"; Axis.side2 <- "right"}
      
      #### Only fit ####
      if(Only.fit == F){FT.line <- geom_line(linetype = "solid", size=0.5, alpha = 0.5)}
      else{
        Geom.point <- NULL
        FT.line <- NULL}
      
      #### Smooth curve ####
      if(Smooth.T == T){
        My_smooth <- geom_smooth(method = "loess", se = Smooth.sd, fullrange = F, level = 0.95, linetype="solid", formula = "y ~ x",
                                 size = .8, show.legend = F, aes(fill = DB), span = Smooth.param[i], alpha = 0.2)}
      else{My_smooth <- NULL}
      
      #### Average lol ####
      if(Mean.models == T){
        Mean.line <- stat_summary(inherit.aes = F, data = Mtot.melt, aes(x = Age, y = variable), geom="line", fun = "mean", color="#770000ff", linewidth = 1.8, alpha = .8, linetype="solid")
      }
      else{Mean.line <- NULL}
      
      #### Plot surf lines + dots ####
      if(is.null(Surf.val.j.i) == F){
        if(is.null(Surf.val.max.age) == T){
          Surf.line <- geom_hline(yintercept = Surf.val.j.i, linetype = "dotdash", size = 0.6, color = "black", alpha = 0.6)
          Surf.dot <- NULL}
        else{
          A = data.frame(x = Limites[1], xend = Surf.val.max.age, yend = Surf.val.j.i, y = Surf.val.j.i)
          Surf.line <- geom_segment(inherit.aes = F, data = A, aes(yend = yend, x = x, xend = xend, y = y), linetype = "dotdash", linewidth = 0.55, color = "black", alpha = 0.55)
          Surf.dot <- geom_point(inherit.aes = F, data = A, aes(x = x, y = y),color = "darkred", alpha = 1, size = Dot.size*1.5, shape = 9)
        }
      }
      else{Surf.dot <- NULL; Surf.line <- NULL}
      
      #### PLOTS ####
      Plot.list[[i]] <- ggplot(data = Mtot.melt, 
                               mapping = aes(x = Age, y = variable, color = DB, shape = Model)) +
        #### Rectangles ####
      scale_fill_manual(values = Rect.color.scale, guide = "none")+
        Rect.clim +
        geom_vline(xintercept = Manual.vlines, col = "grey30", lty = 2, alpha = 0.7)+
        
        #### Plot items ####
      new_scale_fill()+
        scale_shape_manual(values = c("MAT" = 16, "WAPLS" = 1, "RF" = 3, "BRT" = 18)) +
        scale_fill_manual(values = My_color, label = Label.group)+
        scale_color_manual(values = My_color, label = Label.group)+
        Seg.RMSE +
        FT.line + Geom.point + Facet.disp + 
        My_smooth +
        Mean.line + 
        ggtitle(Clim.lab[[j]])+
        guides(colour = guide_legend(nrow = 1))+        # force les legendes a salligner sur une unique ligne
        scale_y_continuous(name = Name.core.plot, breaks = Range.round, position = Axis.side2) +
        scale_x_continuous(name = Title.age, breaks = round(seq(0, Limites[2], by = Select.interv)), expand = c(0.02, 0.02))+
        coord_cartesian(xlim = Limites, ylim = Lim.ano, clip = "off")+
        
        #### Annotations ####
      Repel + Surf.line + Surf.dot +
        
        new_scale_color()+
        scale_color_manual(values = Rect.color.scale, guide = "none", name = NULL, labels = NULL, breaks = NULL, na.translate = FALSE)+
        geom_text(data = yo2, inherit.aes = F, angle = Name.zone.angle,
                  aes(x = (xmax+xmin)/2, y = Inf, label = Title.zone.i, color = Temp.col), size = 3.5, vjust = 1.5, fontface = "bold") +
        ggthemes::geom_rangeframe(data = My_axis_clim, inherit.aes = F, mapping = aes(x = B, y = A), colour = "grey50", sides = Axis.side, size = 0.5) +
        
        #### Theme ####
      theme(
        plot.title =  My_title,
        axis.text.x = element_text(angle = 45, hjust = 1, size = Legende.age.chart[i], colour = "grey55"),
        axis.text.y = element_text(size = 10.5, colour = "grey55"),
        axis.title.x = element_text(size = Legende.age.chart[i], colour = "grey20"),              # taile police du titre de l'axe
        axis.title.y = element_text(size = Legende.Core.name.chart[j], colour = "grey20"),              # taile police du titre de l'axe
        axis.line.y = element_blank(),
        axis.line.x = Xline,
        axis.ticks.x.bottom = Xtick,
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = Legende.position.chart[i],                   # Legendes DB
        legend.justification = c("center"),               # left, top, right, bottom
        legend.direction = "horizontal",
        legend.text.align = 0,
        legend.text = element_text(size = 10, color = "grey20"),
        panel.background = element_blank(),
        panel.spacing = unit(0.08, "cm"),
        panel.grid = element_blank(),
        strip.text.y = element_text(size = 11.5, angle = 180),
        strip.placement = "outside",
        strip.background = element_blank(),
        plot.background = element_blank(),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")
      )
      #### Calcul période climatique ####
      if(Phase.clim == T){
        Continue = NULL
        if(length(unique(Mtot.melt$Model)) == 1 & length(unique(Mtot.melt$DB)) == 1){
          Save.fit <- ggplot_build(Plot.list[[i]])$data[[4]]
          Continue = T}
        if(Mean.models == T){
          Save.fit <- ggplot_build(Plot.list[[i]])$data[[5]]#[,1:6]
          Save.fit <- Save.fit[c("x", "y")]
          Continue = T
        }
        if(is.null(Continue) == T){Continue = F}
        
        if(Continue == T){
          #### New version (loess) ####
          age_grid <- seq(Limites[1], Limites[2], by = 10)
          loess_fit <- loess(variable ~ Age, data = Mtot.melt, span = Smooth.param[i])
          Save.fit <- data.frame(age_grid, predict(loess_fit, newdata = data.frame(Age = age_grid)))
          if(Detrending.fitting == T){Save.fit[[2]] <- resid(lm(Save.fit[[2]] ~ Save.fit[[1]], na.action = na.exclude))}
          Save.fit[[2]] <- 2*(Save.fit[[2]]-min(Save.fit[[2]], na.rm = T))/(max(Save.fit[[2]], na.rm = T)-min(Save.fit[[2]], na.rm = T))-1
          
          #### Clean + export ####
          if(length(unique(Mtot.melt$Model)) == 1 & length(unique(Mtot.melt$DB)) == 1){
            colnames(Save.fit) <- c("Age", paste(Cores[[i]], Pclim[[j]], unique(Mtot.melt$Model), unique(Mtot.melt$DB), sep = "."))
          }
          
          if(Mean.models == T){
            colnames(Save.fit) <- c("Age", paste(Cores[i], Pclim[[j]], "combined_models", sep = "."))
          }
          Merge.mat <- merge(Select.models, Save.fit, by ="Age", all.x = T)
          Select.models <- Merge.mat
          
        }
        
        
      }
      
      #### Verbose fin ####
      if(length(Cores) > 1){
        end[i] <- Sys.time()
        setTxtProgressBar(pb, i)
        time <- round(seconds_to_period(sum(end - init)), 0)
        est <- length(Cores) * (mean(end[end != 0] - init[init != 0])) - time
        remainining <- round(seconds_to_period(est), 0)
        cat(paste(" - Execut. time:", time,
                  " - Estim. time remain.:", remainining), "")}
    }
    if(length(Cores) > 1){close(pb)}
    #### Facet plot clim / cores ####
    Formula.patch <- paste(paste("Plot.list[[", seq(length(Cores)), "]]", sep = ""), collapse = " / ")
    Plot.list.col[[j]] <- eval(parse(text = Formula.patch))
    if(Phase.clim == T){Select.models.tot[[j]] <- Select.models}
  }
  #### Save plots and export ####
  Formula.patch <- paste(paste("Plot.list.col[[", seq(length(Pclim)), "]]", sep = ""), collapse =  " | ")
  Ptot <- eval(parse(text = Formula.patch))
  if(Merge.legends == T){Ptot <- Ptot + plot_layout(guides = "collect")&theme(legend.position = Legend.pos)}
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(Ptot, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", limitsize = F)}
    else{ggsave(Save.plot)}}
  
  if(is.null(Save.Rds) == F & Phase.clim == T){saveRDS(Select.models.tot, Save.Rds)}
  
  if(is.null(Save.plot.Rds) == F){saveRDS(Ptot, Save.plot.Rds)}
  
  return(Ptot)
  
}

Pollen.Ensemble.RMSE <- function(MPal, Plot.x, Keep.clim, Msurf, Mpal.raw = NULL, Mcalib = NULL, Ensemble.WM = T, DES.Ensemble = F, GDGT.plot.merge = F, Variance.method = "quadratic", Temp.col.alpha = 0.1, Cores.lab = NULL,
                                 Limites = NULL, Select.interv = 1000, Lab.age.size = 6, Zone.clim = NULL, Name.zone = NULL, Temp.zone = NULL, Clim.lab = NULL, Dot.size = 6,
                                 Manual.vlines = NULL, Smooth.sd = T, Smooth.T = F, Smooth.param = 0.6, Surf.val = NULL, Name.core = NULL, Surf.val.max.age = NULL, Title.display = T,
                                 Manual.y.val = NULL, Axis.right = F, Anomaly = F, Remove.calib = NULL,
                                 DES.cut.off = NULL, Show.DES.contrib = F, Display.sim = T, DES.Soft.max = 1, DES.distance.metric = "aitchison", DES.KDE = T, DES.parallelize = F,
                                 Add.lim.space = F, Save.plot = NULL, W = NULL, H = NULL, Save.path = NULL, Show.box = F, Plot.Rds = NULL, return.plot = F){
  #### Check if same number of samples in each paleo reconstrcutions ####
  if(length(unique(sapply(MPal, nrow))) > 1){
    print(paste("Be carefull! The ages are not similar from the differentes reconstruction for", Name.core))
    all_ages <- sort(unique(unlist(lapply(MPal, function(df) df$Age))))
    MPal <- lapply(MPal, function(df) {
      template <- data.frame(Age = all_ages)
      out <- full_join(template, df, by = "Age")
      return(out)
    })
  }
  
  #### Weighted variance function ####
  combine_similarity_rmse_models <- function(S, rmses, epsilon = 1e-12) {
    n_pred <- nrow(S)
    datasets <- colnames(S)
    models <- names(rmses)
    
    # Extract dataset name from model names (assumes format "MODEL.DATASET")
    model_datasets <- sub(".*\\.", "", models)
    
    W <- matrix(0, n_pred, length(models))
    
    for (i in 1:n_pred) {
      w_raw <- numeric(length(models))
      
      for (m in seq_along(models)) {
        ds <- model_datasets[m]
        if (!(ds %in% datasets)) next  # skip if dataset not in similarity matrix
        sim_val <- S[i, ds]
        w_raw[m] <- sim_val / (rmses[m]^2 + epsilon)
      }
      
      if (sum(w_raw) == 0) {
        W[i, ] <- rep(1 / length(models), length(models))  # fallback uniform
      } else {
        W[i, ] <- w_raw / sum(w_raw)
      }
    }
    
    colnames(W) <- models
    return(W)
  }
  
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
  
  row_weighted_variance.DES <- function(mat, weights, method = c("standard", "unbiased", "quadratic")) {
    method <- match.arg(method)
    
    if (!all(dim(mat) == dim(weights))) {
      stop("Matrix of predictions and weights must have the same dimensions.")
    }
    
    sapply(seq_len(nrow(mat)), function(i) {
      x <- mat[i, ]
      w <- weights[i, ]
      w_mean <- sum(w * x) / sum(w)
      
      if (method == "standard") {
        return(sum(w * (x - w_mean)^2) / sum(w))
        
      } else if (method == "unbiased") {
        w_sum <- sum(w)
        w_squared_sum <- sum(w^2)
        correction <- w_sum - (w_squared_sum / w_sum)
        return(sum(w * (x - w_mean)^2) / correction)
        
      } else if (method == "quadratic") {
        return(sum(w^2 * (x - w_mean)^2))
      }
    })
  }
  
  Keep.plot.x <- MPal[[1]][Plot.x]
  
  #### Graphical option ####
  Col.pol <- "#BEA33A"
  if(is.null(Limites) == T){Limites = c(min(MPal[[1]][Plot.x], na.rm = T), max(MPal[[1]][Plot.x], na.rm = T))}
  
  if(Show.box == T){
    My_border <- element_rect(fill = NA, colour = "grey30", linewidth = .3)
    Axis.x <- element_blank(); Axis.y <- element_blank()
  }
  else{
    My_border <- element_blank()
    Axis.x <- element_line(); Axis.y <- element_line()
  }
  
  if(GDGT.plot.merge == T){
    Title.age = NULL
    Axis.x <- element_blank(); Ticks.x <- element_blank(); Axis.lab.x <- element_blank()
    if(Add.lim.space == T){
      Limites[2] <- (Limites[2] + Limites[2]*0.23)}
    Legende.Core.name.chart <- c(17, rep(0, (length(Keep.clim)-1)))
  }
  else{
    Title.age = "Age (yr cal BP)"
    Ticks.x <- element_line(); Axis.lab.x <- element_text(angle = 45, hjust = 1, size = Lab.age.size, colour = "grey55")
    Legende.Core.name.chart <- c(16, rep(0, (length(Keep.clim)-1)))
  }
  
  if(Title.display == F){My_title <- element_blank()}
  else{My_title <- element_text(size = 17, hjust = 0.5)}
  
  Plot.list.col <- list()
  
  My_color_models <- c("WASTDB" = "#e2a064ff", "WAST" = "#e2a064ff", 
                       "NMSDB" = "#F3A481", "MDB" = "#963326", "TUDB" = "#0094AF", "TUSDB" = "#0094AF",
                       "ST" = "#91C4DD", "STDB" = "#91C4DD", 
                       "ACADB" = "#963326", "CAUCDB" = "#0F3361",
                       "COSTDB" = "#f3c768ff", "COST" = "#f3c768ff",
                       "MEDTEMP" = "#F3A481", "TEMPSCAND" = "#91C4DD", 
                       "EAPDB" = "#0F3361", "TAIGDB" = "#32156eff")
  
  #### Climat zone settings ####
  Rect.color.scale <- c("W" = "#E76D51",
                        "C" = "#75AADB",
                        "G" = "grey40",
                        "D" = "#c67f05",
                        "Wt" = "#004266")
  Rect.color.scale <- Rect.color.scale[unique(Temp.zone)]
  
  if(is.null(Zone.clim) ==F){
    yo = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                    xmax = Zone.clim[seq(2,length(Zone.clim), by=2)], 
                    Temp.col = Temp.zone)
    Rect.clim <- geom_rect(data = yo, inherit.aes = F,
                           mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf, fill = Temp.col),
                           alpha = Temp.col.alpha, color = "grey", size = 0.3, linetype = 2, na.rm = T)} 
  else{Rect.clim <- NULL}
  if(is.null(Name.zone) == F){yo2 = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                                               xmax = Zone.clim[seq(2,length(Zone.clim), by=2)],
                                               Temp.col = Temp.zone,
                                               Title.zone = Name.zone,
                                               Model = Model[length(Model)])}
  else{yo2 = data.frame(xmin = 0, xmax = 0, Temp.col = "", Title.zone = "")}
  
  #### Case only some Surface Values ####
  if(is.null(Surf.val) == F){
    A <- setNames(data.frame(matrix(ncol = length(setdiff(Keep.clim, names(Surf.val))), nrow = nrow(Surf.val))), setdiff(Keep.clim, names(Surf.val))) 
    row.names(A) <- row.names(Surf.val)
    A[is.na(A)] <- NA
    A <- cbind(Surf.val, A)
    Surf.val <- A[,sort(names(A))]}
  
  #### Smooth curve ####
  if(Smooth.T == T){
    My_smooth <- geom_smooth(method = "loess", se = Smooth.sd, fullrange = F, level = 0.95, linetype="solid", formula = "y ~ x",
                             size = .8, show.legend = F, fill = Col.pol, color = Col.pol, span = Smooth.param, alpha = 0.2)}
  else{My_smooth <- NULL}
  #### Dynamic Ensemble Selection based on Aitchison distance ####
  if(DES.Ensemble == T){
    if(is.null(Mcalib) == T | is.null(Mpal.raw) == T){
      print("**** Dynamic Ensemble Selection impossible without 'Mcalib' and 'Mpal.raw'. ****")
      Msim <- NULL; DES.Ensemble = F}
    else{
      Keep.calib <- unique(gsub(".*\\.", "", names(MPal)))
      Keep.calib2 <- unique(gsub(".*\\.", "", names(Msurf)))
      Keep.calib <- intersect(Keep.calib, Keep.calib2)
      
      if(length(names(Mcalib)[!names(Mcalib) %in% Keep.calib]) > 0){
        
        print(paste("**** The following dalibrations were removed from the DES, since they were not applied on ", Name.core, ". ****", sep = ""))
        print(names(Mcalib)[!names(Mcalib) %in% Keep.calib])
      }
      
      
      if(all(Keep.calib %in% names(Mcalib)) == F){
        print(paste("**** Dynamic Ensemble Selection impossible without ", Keep.calib[!Keep.calib %in% names(Mcalib)], " given in 'Mcalib'. ****", sep = ""))
        Msim <- NULL
        DES.Ensemble = F
      }
      else{
        Mcalib <- Mcalib[names(Mcalib) %in% Keep.calib]
        
        if(unique(lapply(MPal, nrow)) != nrow(Mpal.raw)){
          print("**** The number of sample in 'MPal' and 'Mpal.raw are different ! ****")
          Msim <- NULL; Msim.m <- NULL; Show.DES.contrib = F; Display.sim = F; DES.Ensemble = F
        }
        else{
          Msim <- DES.weights(
            models_cal = Mcalib, Xpred = Mpal.raw, 
            K_cal = min(sapply(Mcalib, nrow)), K_pred = signif(nrow(Mpal.raw)/10, digits = 0),
            Cut.off = DES.cut.off, Soft.max = DES.Soft.max, causal = F, n_cores = 2, 
            metric = DES.distance.metric, Local.subsampling = DES.KDE, parallelize = DES.parallelize, )
          
          if(Show.DES.contrib == T & Display.sim == T){
            Msim.m <- cbind(Keep.plot.x, Msim)
            Msim.m <- melt(Msim.m, id = Plot.x)
            names(Msim.m)[1] <- "Plot.x"
          }
        }
      }
    }
  }
  else{Msim <- NULL; Msim.m <- NULL}
  
  #### Loop on climate param ####
  Res <- list()
  Res.contrib <- list()
  for(j in Keep.clim){
    List.predict <- list(); RMSE <- c(NA)
    for(i in 1:length(Msurf)){
      Name.model.full <- paste(j, names(Msurf)[i], sep = ".")
      
      Method.i <- names(Msurf)[i]
      if(is.null(Remove.calib) == F){
        print(any(Remove.calib %in% Name.model.full))
        # if(Remove.calib == Name.model.full){next}
        if(any(Remove.calib %in% Name.model.full)){next}
      }
      
      if(j %in% names(MPal[[Method.i]]) == F){next}
      
      List.predict[[i]] <- setNames(MPal[[Method.i]][j], Method.i)
      RMSE.i <- Msurf[[Method.i]]$Best.Param
      RMSE.i <- RMSE.i[row.names(RMSE.i) == j, "RMSE.choice"]
      RMSE[i] <- RMSE.i
    }
    List.predict <- Filter(Negate(is.null), List.predict)
    RMSE <- Filter(Negate(is.na), RMSE)
    List.predict <- do.call(cbind, List.predict)
    
    #### Weighted Mean RMSE ####
    if(Ensemble.WM == T){
      #### DES ####
      if(DES.Ensemble == T){
        names(RMSE) <- names(List.predict)
        wg <- combine_similarity_rmse_models(Msim, RMSE)
        
        Ensemble_mean <- rowSums(List.predict * wg)
        prediction_matrix <- as.matrix(List.predict)
        weighted_variance <- row_weighted_variance.DES(prediction_matrix, wg, method = Variance.method)
        
        Res[[paste(j, "FT_Ensemble", sep = "_")]]    <- Ensemble_mean
        Res[[paste(j, "FT_Ensemble_MaxI", sep = "_")]] <- Ensemble_mean + 1.96 * sqrt(weighted_variance)
        Res[[paste(j, "FT_Ensemble_MinI", sep = "_")]] <- Ensemble_mean - 1.96 * sqrt(weighted_variance)
        
        if(Show.DES.contrib == T & Display.sim == F){
          colnames(wg) <- paste(j, colnames(wg), sep = "_")
          Res.contrib[[j]] <- wg
        }
      }
      
      #### Only RMSE ####
      if(DES.Ensemble == F){
        wg <- 1/RMSE
        wg <- wg/sum(wg, na.rm = T)
        Res[[paste(j, "FT_Ensemble", sep = "_")]] <- as.numeric(as.matrix(List.predict) %*% wg)
        prediction_matrix <- as.matrix(List.predict)
        weighted_variance <-row_weighted_variance(prediction_matrix, wg, method = Variance.method)
        Res[[paste(j, "FT_Ensemble_MaxI", sep = "_")]] <- Res[[paste(j, "FT_Ensemble", sep = "_")]] + 1.96*sqrt(weighted_variance)
        Res[[paste(j, "FT_Ensemble_MinI", sep = "_")]] <- Res[[paste(j, "FT_Ensemble", sep = "_")]] - 1.96*sqrt(weighted_variance)
      }
    }
  }
  Res <- cbind(Keep.plot.x, do.call(cbind, Res))
  
  #### Prep melted data ####
  
  melt <- melt(Res, names(Res)[c(1, grep("MinI", names(Res)), grep("MaxI", names(Res)))])
  names(melt)[c(ncol(melt)-1, ncol(melt))] <- c("Param.clim", "Mean_val")
  
  melt <- melt(melt, c("Param.clim", "Mean_val", names(melt)[c(1, grep("MaxI", names(melt)))]))
  names(melt)[c(ncol(melt)-1, ncol(melt))] <- c("MinI", "Min_val")
  
  melt <- melt(melt, c("Param.clim", "Mean_val", "MinI", "Min_val", Plot.x))
  names(melt)[c(ncol(melt)-1, ncol(melt))] <- c("MaxI", "Max_val")
  
  names(melt)[names(melt) == Plot.x] <- "Plot.x"
  melt$Param.clim <- gsub("_FT_Ensemble", "", melt$Param.clim)
  melt$MaxI <- gsub("_FT_Ensemble_MaxI", "", melt$MaxI)
  melt$MinI <- gsub("_FT_Ensemble_MinI", "", melt$MinI)
  
  melt <- melt[melt$Param.clim == melt$MinI & melt$Param.clim == melt$MaxI,]
  melt$Param.clim <- factor(melt$Param.clim, levels = Keep.clim, ordered = T)
  melt <- melt[!duplicated(melt),]
  
  #### Melt contrib ####
  if(Show.DES.contrib == T & Display.sim == F){
    Res.contrib <- cbind(Keep.plot.x, do.call(cbind, Res.contrib))
    Res.contrib <- melt(Res.contrib, Plot.x)
    
    print(Res.contrib)
    # A finir !
    
  }
  
  #### Plot ####
  for(i in 1:length(Keep.clim)){
    melt.i <- melt[melt$Param.clim == Keep.clim[i],]
    if(is.null(Surf.val) == F){Surf.val.j <- subset(Surf.val, select = c(Keep.clim[i]))}
    
    #### Surf val ####
    if(is.null(Surf.val) == F){
      Surf.val.j.i <- Surf.val.j[Name.core,]
      if(is.na(Surf.val.j.i) == T){Surf.val.j.i = NULL}
    }
    
    if(is.null(Surf.val) == F){
      if(is.na(Surf.val.j.i) == F){
        if(is.null(Surf.val.max.age) == T){
          Surf.line <- geom_hline(yintercept = Surf.val.j.i, linetype = "dotdash", linewidth = 0.55, color = "black", alpha = 0.55)
          Surf.dot <- NULL}
        else{
          A = data.frame(x = Limites[1], xend = Surf.val.max.age, yend = Surf.val.j.i, y = Surf.val.j.i)
          Surf.line <- geom_segment(inherit.aes = F, data = A, aes(yend = yend, x = x, xend = xend, y = y), linetype = "dotdash", linewidth = 0.55, color = "black", alpha = 0.55)
          Surf.dot <- geom_point(inherit.aes = F, data = A, aes(x = x, y = y),color = "darkred", alpha = 1, size = Dot.size*1.5, shape = 9)
        }
        
      }
      if(is.na(Surf.val.j.i) == T | is.null(Surf.val.j.i) == T){Surf.line <- NULL; Surf.dot <- NULL}}
    else{Surf.line <- NULL; Surf.dot <- NULL}
    
    #### Breaks ####
    if(is.null(Manual.y.val) == T){
      if(is.null(Surf.val) == F){
        Kmax <- max(melt.i$variable, Surf.val.j.i, na.rm = T)
        Kmin <- min(melt.i$variable, Surf.val.j.i, na.rm = T)}
      else{
        Kmax <- max(melt.i$variable, na.rm = T)
        Kmin <- min(melt.i$variable, na.rm = T)
      }
      
      if(Kmax >= 40){
        Kmin <- round(Kmin, digits = -1)
        Kmin = floor(Kmin / 50) * 50
        Kmax <- round(Kmax, digits = -1)
        Kmax = ceiling(Kmax / 50) * 50
        Kmid <- Kmin + (Kmax - Kmin)*1/3
        Kmid = ceiling(Kmid / 50) * 50
        Kmid2 <- Kmin + (Kmax - Kmin)*2/3
        Kmid2 = ceiling(Kmid2 / 50) * 50
      }
      if(Anomaly == T & is.null(Surf.val) == F){
        Kmid <- Kmin + (Kmax - Kmin)*1/3
        Kmid2 <- Kmin + (Kmax - Kmin)*2/3
      }
      else{
        Kmin <- round(Kmin, digits = 1)
        Kmin = floor(Kmin / 0.5) * 0.5
        Kmax <- round(Kmax, digits = 1)
        Kmax = ceiling(Kmax / 0.5) * 0.5
        Kmid <- Kmin + (Kmax - Kmin)*1/3
        Kmid = ceiling(Kmid / 0.5) * 0.5
        Kmid2 <- Kmin + (Kmax - Kmin)*2/3
        Kmid2 = ceiling(Kmid2 / 0.5) * 0.5
      }
      Range.round <- c(Kmin, Kmid, Kmid2, Kmax)
      if(Kmax >= 40){
        Range.round <- round(Range.round, digits = -1)
      }
      if(Anomaly == T & is.null(Surf.val) == F){
        Range.round <- round(Range.round, digits = 2)
      }
      else{
        Range.round <- round(Range.round, digits = 1)
      }
      
    }
    else{
      Clim.break <- gsub("\\..*", "", names(Manual.y.val))
      Core.break <- gsub(".*\\.", "", names(Manual.y.val))
      Match.clim <- which(Clim.break %in% Keep.clim[i])
      Match.tot <- unique(Match.clim)
      Range.round <- Manual.y.val[[Match.tot]]
    }
    
    My_axis_clim <- data.frame(A = Range.round[c(1,length(Range.round))], B = Limites)
    if(Axis.right == F){Axis.side <- "l"; Axis.side2 <- "left"}
    else{Axis.side <- "r"; Axis.side2 <- "right"}
    
    if(is.null(Cores.lab) == F){
      if(is.null(Manual.y.val) == T){Scale_y <- scale_y_continuous(name = Cores.lab,  breaks = scales::breaks_extended(n = 4))}
      else{Scale_y <- scale_y_continuous(name = Cores.lab, breaks = Range.round, position = Axis.side2)}
    }
    else{Scale_y <- NULL}
    
    #### Plot ####
    p <- ggplot(melt.i, aes(x = Plot.x, y = Mean_val))+
      #### Rectangles ####
    scale_fill_manual(values = Rect.color.scale, guide = "none")+
      Rect.clim +
      geom_vline(xintercept = Manual.vlines, col = "grey30", lty = 2, alpha = 0.7)+
      new_scale_fill()+
      
      #### Items ####
    facet_wrap(vars(Param.clim), scale = "free")+
      geom_ribbon(aes(ymin = Min_val, ymax = Max_val), fill = Col.pol, color = NA, size = 0.15, linetype = "dashed", alpha = 0.1)+
      
      xlab(Plot.x)+
      geom_point(size = Dot.size, color = Col.pol, shape = 1)+
      geom_line(color = Col.pol, linetype = "solid", linewidth = 0.5, alpha = 0.6)+
      scale_x_continuous(name = Title.age, breaks = round(seq(0, Limites[2], by = Select.interv)), expand = c(0.02, 0.02))+
      Scale_y +
      My_smooth +
      ggtitle(Clim.lab[[i]])+
      guides(colour = guide_legend(nrow = 1))+        # force les legendes a salligner sur une unique ligne
      #### Annotations ####
    new_scale_color()+
      scale_color_manual(values = Rect.color.scale, guide = "none", name = NULL, labels = NULL, breaks = NULL, na.translate = FALSE)+
      # geom_text(data = yo2, inherit.aes = F,
      # aes(x = (xmax+xmin)/2, y = Inf, #label = Title.zone.i, 
      # color = Temp.col), size = 3.5, vjust = 1.5, fontface = "bold") +
      ggthemes::geom_rangeframe(data = My_axis_clim, inherit.aes = F, mapping = aes(x = B, y = A), colour = "grey50", sides = Axis.side, size = 0.5) +
      Surf.line + Surf.dot +
      #### Theme ####
    theme(plot.background = element_blank(), axis.ticks.x = Ticks.x,
          plot.title = My_title, panel.grid = element_blank(),
          element_text(size = 11, colour = "grey20"),
          axis.title.y = element_text(size = Legende.Core.name.chart[i], colour = "grey20"),              # taile police du titre de l'axe
          strip.text = element_blank(),
          axis.line.x = Axis.x, axis.line.y = Axis.y,
          panel.background = My_border, strip.background = element_blank(),
          axis.text.x = Axis.lab.x,
          plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
          
          legend.position = "none")
    
    
    #### Plot contribution DES ####
    if(Show.DES.contrib == T & DES.Ensemble == T){
      p2 <- ggplot(Msim.m, aes(x = Plot.x, y = value, color = variable, fill = variable)) +
        geom_col(position = "stack") +
        scale_color_manual(values = My_color_models)+
        scale_fill_manual(values = My_color_models)+
        scale_x_continuous(name = Title.age, breaks = round(seq(0, Limites[2], by = Select.interv)), expand = c(0.02, 0.02))+
        theme(plot.background = element_blank(), 
              plot.title = element_blank(), panel.grid = element_blank(),
              element_text(size = 11, colour = "grey20"),
              strip.text = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
              panel.background = My_border, strip.background = element_blank(),
              plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
              legend.position = "top", legend.direction = "horizontal", legend.title = element_blank(),
        )
      
      p <- p2 / p + plot_layout(heights = c(25,75))
    }
    
    Plot.list.col[[i]] <- p
  }
  #### Facet plot clim / cores ####
  Formula.patch <- paste(paste("Plot.list.col[[", seq(length(Keep.clim)), "]]", sep = ""), collapse =  " | ")
  p <- eval(parse(text = Formula.patch))
  
  #### Export plots #### 
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(p, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", limitsize = F)}
    else{ggsave(Save.plot)}}
  
  if(is.null(Plot.Rds) == F){saveRDS(p, Plot.Rds)}
  if(is.null(Save.path) == F){saveRDS(Res, Save.path)}
  if(return.plot == T){return(p)}
  else{return(Res)}
}

Timeline.build <- function(Datatime = NULL, Limites = NULL, Axis.pos.y = 0, Select.interv = 1000, Label.size = 4.2, 
                           Add.AC.BC = F, 
                           Leg.pos = "none", Save.plot = NULL, W = NULL, H = NULL, Save.plot.Rds = NULL) {
  #### Clean data ####
  if(is.null(Limites) == F){
    Datatime$Start[Datatime$Start > Limites[2]] <- Limites[2]
    My_breaks <- seq(0, max(Datatime$Start), Select.interv)
    Val.AD <- 2000 - My_breaks
  }
  
  Datapoint <- Datatime[Datatime$Span == "Event",]
  Datatime <- Datatime[Datatime$Span == "Period",]
  
  #### Graphical param ####
  if(Add.AC.BC == T){
    My_ticks <- element_line()
    My_axis_text <- element_text(size = 8)
    Axis.x <- scale_x_continuous(name = "Time (year BC/AD)", breaks = My_breaks, expand = c(0.02, 0.02), labels = Val.AD)
    My_axis_tit <- element_text()
  }
  else{
    Axis.x <- scale_x_continuous(expand = c(0.02, 0.02))
    My_ticks <- element_blank(); My_axis_text <- element_blank(); My_axis_tit <- element_blank()
  }
  print(Datapoint)
  #### Plot ####
  time_plot <- Datatime %>% 
    ggplot(aes(x = Start, y = Type)) +
    geom_segment(data = Datapoint, mapping =  aes(x = Start, xend = Start, y = Type), yend = Axis.pos.y, color='grey30', size = Label.size/8, linetype = "dashed", alpha = .7, show.legend = T) +
    geom_segment(y = Axis.pos.y, yend = Axis.pos.y, x = Limites[1], xend = max(Datatime$Start), color='grey30', size = Label.size/8, alpha = .7, show.legend = F) +
    
    geom_segment(aes(xend = End, yend = Type, color =  Carac, size = Seg_size)) +
    
    geom_text(aes(x = (Start + End)/2, label = Lab, vjust = Vjust, color = Carac_lab, angle = Lab_ang, hjust = Hjust), size = Label.size, show.legend = F) +
    
    scale_color_manual(values = c("Wet" = "#004266", "Urban" = "#D62828", "Warm" = "#E76D51", "Cold" = "#75AADB", "Dry" = "#c67f05", "Nomad" = "#03875B", "nul" = "white"))+ #FCB322
    scale_size_continuous(range = c(Label.size/2.25,Label.size*1.7))+
    geom_point(data = Datapoint, mapping = aes(x = Start, y = Type, color = Carac), size = Label.size/2, show.legend = F)+
    geom_text(data = Datapoint, mapping = aes(x = Start, label = Lab, hjust = ifelse(Vjust > 0, 1.1, -.1), color = Carac, angle = Lab_ang),
              vjust = -.4, size = Label.size, show.legend = F) +
    ylim(c(0,6))+
    coord_cartesian(clip = 'off') +
    
    Axis.x +
    guides(fill = "none", size = "none", color = guide_legend(nrow = 1)) +
    theme(legend.position = Leg.pos,
          legend.title = element_blank(),
          legend.direction = "horizontal",
          axis.title.x = My_axis_tit,
          axis.title.y=element_blank(), panel.grid = element_blank(),
          axis.text.x = My_axis_text,
          axis.text.y = element_blank(),
          axis.ticks.y =element_blank(),
          axis.ticks.x = My_ticks,
          axis.line.x =element_blank(),
          # strip.clip = "off", plot.margin = unit(c(0,0,0,0), 'cm'),
          panel.background = element_blank(),
          plot.background = element_blank(),
          text = element_text(size = 20))
  
  #### Save plot and export ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(time_plot, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.plot)}}
  
  if(is.null(Save.plot.Rds) == F){saveRDS(time_plot, Save.plot.Rds)}
  return(time_plot)
}
