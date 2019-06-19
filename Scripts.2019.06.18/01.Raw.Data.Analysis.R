######################################################################################################
# 01.Raw.Data.Analysis.R
######################################################################################################
try.dev.off(); # rm(list=ls(all.names = TRUE));
# source("~/GitHub/Kinetochore/Scripts/01.Raw.Data.Analysis.R")


# Setup ------------------------

InputDir = p0(InputDirOrig, "2017.08.23/" )
setup_MarkdownReports(OutDir = paste0(OutDirOrig,"01.Raw.Data.Normalisation"), scriptname = "01.Raw.Data.Analysis.R")
LS_prots = list.files(path = InputDir, pattern = "^[^~]*.xlsx$", all.files	=F) # ignore Find files starting with ~$ (MicroSoft Word Temp files)

TheProteins = stringr::str_split_fixed(LS_prots,pattern = '.xlsx$', n=2)[,1]


lapply(TheProteins, data.frame)
# Parameters ------------------------
ReadIn = T

# Read in ------------------------
try.dev.off()
if (ReadIn) {
  LS_excels = list.fromNames(LS_prots)
  for (i in 1:length(LS_prots) ) {
    i
    path =kollapse(InputDir,LS_prots[i])
    file.exists(path)
    LS_excels[[i]] = read.simple.xls(pfn = path)
  } #for
  LS_excels.bac = LS_excels
} else {
  LS_excels = LS_excels.bac
} # if (ReadIn) 
names(LS_excels)

# Prepare Clean table ------------------------
CENPC_norm_ls = list.fromNames(SHEETZ); MATnames= NULL
for (i in seq(from = 3,to = 21, by = 3) ) {
  MATnames[i-2] = TheProteins[i/3]
  MATnames[i-1] = p0("TUB.", TheProteins[i/3])
  MATnames[i] = p0("BATCH.", TheProteins[i/3])
} #for
MATnames
CENPC_norm_ls = lapply(CENPC_norm_ls, matrix.fromNames, rowname_vec = 1:200, colname_vec = MATnames)
CENPC_norm_ls = lapply(CENPC_norm_ls, as.data.frame)

# ------------------------------------------------------------------------------------------

# "Fig.S3D-E ------------------------------------------------------------------------------------------
# i=1
pdfA4plot_on(pname = "Fig.S1D-E.RawIntesisites.signal.2.BG", cols = 4)
for (i in 1:length(LS_prots) ) {
  p = LS_prots[i]
  protname = strsplit(p,"\\.")[[1]][1]
  dat = LS_excels[[p]]
  stopifnot(SHEETZ == names(dat))
  
  for (j in 1:length(SHEETZ) ) {
    e=SHEETZ[j]
    sht = dat[[e]]
    sign2noise = T
    if (sign2noise) {
      sht$"s2n.TUB.log2" = log2(sht[ ,"TUB"] / sht[ ,"BG"])
      sht$"s2n.prot.log2" = log2(sht[ ,3] / sht[ ,"BG.1"])
      sht$"s2n.CENPC.log2" = log2(sht[ ,"CENP.C"] / sht[ ,"BG.2"])
    } #if
    COLUMNZ = colnames(sht)
    PlottingOrder = c(1,2,8,11,  3,4,9,12,    5,6,10,13 )
    COL ='s2n.prot.log2'
    ls_raw_per_exp = split(sht[,COL] , f = sht[,"BATCH"])
    names(ls_raw_per_exp) = p0('e',names(ls_raw_per_exp))
    
    ylabel = "log2(Signal/Background)" 
    wstripchart(ls_raw_per_exp, savefile = F, pch = 18, pchcex = 1, ylb = ylabel, main = paste(protname, e))
    
    # CENPC_norm ------------------------------------------------------
    CENPC_norm = T
    if (CENPC_norm) {
      # View(LS_excels[[p]][[e]])
      LS_excels[[p]][[e]]$"CENPC_norm.TUB" = sht$"CENPC_norm.TUB" = (sht[ ,"TUB.BG"] / sht[ ,"CENP.C.BG"])
      LS_excels[[p]][[e]]$"CENPC_norm.prot2" = sht$"CENPC_norm.prot2" = (sht[ ,9] / sht[ ,"CENP.C.BG"]); COLUMNZ[9]
      
      ls_CENPC_norm.TUB = split(sht[,"CENPC_norm.TUB"] , f = sht[,"BATCH"])
      ls_CENPC_norm.prot2 = split(sht[,"CENPC_norm.prot2"] , f = sht[,"BATCH"])
    } #if CENPC_norm
  } #for
  plot.new()
} #for
pdfA4plot_off()


# Fig.S3F ------------------------
for (j in 1:length(SHEETZ) ) {
  e=SHEETZ[j]
  pname = p0("Fig.S3F.CENPC.correlation.", e)
  pdfA4plot_on(pname = pname, rows = 5)
  for (i in 1:length(LS_prots) ) {
    p = LS_prots[i]
    dat = LS_excels[[p]]
    SHEETZ = names(dat)
    sht = dat[[e]]

    ccc = rich.colors.vec(as.factor.numeric(sht$BATCH))
    ProtX = colnames(sht)[3]
    Prot_and_BG= sht[ ,c(ProtX, "CENP.C.BG")]
    colnames(Prot_and_BG)= c(p0(ProtX, "-BG [EXP]"), "CENP-BG [EXP]")

    plot(Prot_and_BG, bg = ccc, pch = 21)
  } #for EXCEL FILE / PROTEIN OF INTEREST
  pdfA4plot_off()
} #for SHEET


# ---------------------------------------------------------------------------------------------------------

writeOutNormData=T
if (writeOutNormData) {
  # compile list of tables: FULL, NULL, EXP each 3 columns * each protein
  for (i in 1:length(LS_prots) ) {
    p = LS_prots[i]
    prot = strsplit(p,"\\.")[[1]][1]
    dat = LS_excels[[p]]
    SHEETZ = names(dat)
    for (j in 1:length(SHEETZ) ) {
      e=SHEETZ[j]; e
      NrMeasm=NROW(LS_excels[[p]][[e]])
      CENPC_norm_ls[[e]][ , prot               ][1:NrMeasm] = LS_excels[[p]][[e]]$"CENPC_norm.prot2"
      CENPC_norm_ls[[e]][ , p0("TUB.", prot)   ][1:NrMeasm] = LS_excels[[p]][[e]]$"CENPC_norm.TUB"
      CENPC_norm_ls[[e]][ , p0("BATCH.", prot) ][1:NrMeasm] = LS_excels[[p]][[e]]$"BATCH"
    }
  } # for length(LS_prots)
  
  for (j in 1:length(SHEETZ) ) {
    e=SHEETZ[j]
    fname = p0(InputDir, e, ".tsv"); fname
    write.simple.tsv(CENPC_norm_ls[[e]], ManualName = fname)
  }
} #if


# Normalization calc from raw data ----------------------------------------------------------------------------------------------------

Final.Data = list.fromNames(LS_prots)
for (i in 1:length(LS_prots) ) {
  p = LS_prots[i]
  dat = LS_excels[[p]]
  SHEETZ = names(dat)
  
  nData = dat$`NULL`
  fData = dat$FULL
  eData = dat$EXP
  
  # 2 variables for FULL (tub) and NULL (prot) normalisation
  fData$"CENPC_norm.TUB" = (fData[ ,"TUB.BG"] / fData[ ,"CENP.C.BG"])
  nData$"CENPC_norm.prot2" = (nData[ ,9] / nData[ ,"CENP.C.BG"])
  
  # Variables to be normalised
  eData$"CENPC_norm.TUB" = (eData[ ,"TUB.BG"] / eData[ ,"CENP.C.BG"])
  eData$"CENPC_norm.prot2" = (eData[ ,9] / eData[ ,"CENP.C.BG"])
  
  med.TUB.per.Exp.FULL = unlapply(split(fData$"CENPC_norm.TUB", f = fData$BATCH), median)
  med.Prot.per.Exp.NULL = unlapply(split(nData$"CENPC_norm.prot2", f = nData$BATCH), median)
  
  median.per.Exp.TUB.Full = med.TUB.per.Exp.FULL[eData$BATCH]
  median.per.Exp.Prot.Null = med.Prot.per.Exp.NULL[eData$BATCH]
  
  TUB.FullNormalised = eData$"CENPC_norm.TUB" / median.per.Exp.TUB.Full
  Prot.NullNormalised = eData$"CENPC_norm.prot2" / median.per.Exp.Prot.Null
  Final.Data[[i]] = Ready.EXP = cbind(Prot.NullNormalised, TUB.FullNormalised )
  BATCHES = rownames(Ready.EXP) = eData$BATCH
  Ready.EXP.per.EXP = split.data.frame(Ready.EXP, f = BATCHES)
  
} #for EXCEL FILE / PROTEIN OF INTEREST

CombineFinal2DF = T
if (CombineFinal2DF) {
  Final.Data.DF =list2df(Final.Data)
  # new names
  NameParts = (stringr::str_split_fixed(colnames(Final.Data.DF), pattern = "\\.", n = 4)[, c(1,3)])
  NewNames = paste(NameParts[,1], NameParts[,2], sep = ".")
  colnames(Final.Data.DF) = NewNames
} #if


