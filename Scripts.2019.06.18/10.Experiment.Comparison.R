################################################################
# 10.Experiment.Comparison
################################################################
# source("~/GitHub/Kinetochore/Scripts/10.Experiment.Comparison.R")


# Setup ------------------------
setup_MarkdownReports(OutDir = p0(OutDirOrig, "10.Experiment.Comparison/"), scriptname = "10.Experiment.Comparison.R")

# Parameters ------------------------


# Data ------------------------
InDir = InputDir
Exp = read.simple.tsv(p0(InDir, "EXP.tsv" ), NaReplace = F)
Full = read.simple.tsv(p0(InDir, "FULL.tsv"), NaReplace = F )
Null = read.simple.tsv(p0(InDir, "NULL.tsv"), NaReplace = F )

any(colnames(Null) != colnames(Full))
woBatchCol = grepv(pattern = "^BATCH" , colnames(Null), invert = T)

Exp = Exp[, woBatchCol]
Full = Full[, woBatchCol]
Null = Null[, woBatchCol]

# Meta data ------------------------
Measurements = colnames(Null)

Proteins = grep("^TUB", Measurements, invert = T, value = T); Proteins
K.fibers = grep("^TUB", Measurements, invert = F, value = T); K.fibers


# --------------------------------------------------------------------------------
# Proteins
Proteins_FULL = as.list.data.frame(Full[,Proteins])
Proteins_NULL = as.list.data.frame(Null[,Proteins])

Medianz_N = unlapply(Proteins_NULL, median, na.rm = T)
Medianz_F = unlapply(Proteins_FULL, median, na.rm = T)
ProteinDecreaseOnFullAttachment = (100*Medianz_F / Medianz_N)

# Plots --------------------------------------------------------------------------------

SD_F = unlapply(Proteins_FULL, sd, na.rm = T)
Fig.S3C.Decrease.and.SD =  cbind(
  "Standard Deviation (%)" = 100*SD_F,
  "Relative Decresease % ( NULL/FULL)" = ProteinDecreaseOnFullAttachment)
stopif(sum(names(SD_F) != names(ProteinDecreaseOnFullAttachment)), message = "The (order) of the names does not match.")

# ccc= rich.colors.vec((names(SD_F)))
ccc = c(ZW10 = "#F6F906FF", Spindly = "#13F24AFF", pMELT = "#FFBA00FF", 
        MAD2 = "#00004AFF", MAD1 = "#0010FFFF", BUBR1 = "#00A4FFFF", 
        BUB1 = "#FF3300FF")
# color_check(ccc)

wplot(Fig.S3C.Decrease.and.SD, cex=2, type = "n")
abline(h=seq(0,50,10), lty=3, col="grey"); abline(v=seq(0,18,2), lty=3, col="grey")
points(Fig.S3C.Decrease.and.SD, pch=23, bg=ccc[rownames(Fig.S3C.Decrease.and.SD)], cex=2)
wlegend(sortbyitsnames(ccc), 1)

# --------------------------------------------------------------------------------




