######################################################################
# 05.Fig3.Tests.Normality.Variance.R
######################################################################
# source ('~/GitHub/Kinetochore/Scripts/05.Fig3.Tests.Normality.Variance.R')
# rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
irequire(car)
irequire(pheatmap)
# require(data.table)
library(reshape2)

# Setup ------------------------

setup_MarkdownReports(OutDir = p0(OutDirOrig,"05.Fig3.Distributions.and.variance"),
                      scriptname = "05.Fig3.Tests.Normality.Variance.R", title = "Tests.Normality.Variance")

# InDir = p0(InputDirOrig, "/2017.06.06/")
usenew =T

InDir = "~/GitHub/Kinetochore/Data/2017.08.23/"
Full = read.simple.tsv(p0(InDir, "FULL.tsv"), NaReplace =  F)
Null = read.simple.tsv(p0(InDir, "NULL.tsv"), NaReplace =  F)

any(colnames(Full) != colnames(Null))
woBatchCol = grepv(pattern = "^BATCH", colnames(Full), invert = T)

Full = Full[, woBatchCol]
Null = Null[, woBatchCol]



NrProteins = NCOL(Null)/2

n=4
set.seed(21768134)
ccc = sample(richColors(2*NrProteins+n)[-(1:n)] )

# Data normalisation ------------------------
llprint("## Data normalisation for Fig3")
llprint("> As of Thu Jul 13 21:17:46 2017")
llprint("Both Full and Null (-attachment) datasets are normalised to:")
llprint("- their local background (in the same channel, done before)")
llprint("- CENP-C signal at the same locus")
llprint("- Division by the **median** of the Null datased")
llprint("    - Null dataset will be centered around 1")
llprint("    - Full dataset will represent relative concentraition")

Null.med = colMedians(Null);l(Null.med)

Full.divbyNull = colDivide(Full , Null.med)
Null.norm.self = colDivide(Null, Null.med)

Full.divbyNull.Prot = get.oddoreven(Full.divbyNull, odd = usenew)
Null.divbyNull.Prot = get.oddoreven(Null.norm.self, odd = usenew)
Protz.NullNorm = intermingle.cbind(Null.divbyNull.Prot, Full.divbyNull.Prot )

# Basic statistics ------------------------------------------------------------------------------------------
Stat.prots = iround(cbind(
  "var" =  apply(Full.divbyNull.Prot, 2, var, na.rm=T),
  "std" =  apply(Full.divbyNull.Prot, 2, sd, na.rm=T),
  "cv" =   apply(Full.divbyNull.Prot, 2, cv, na.rm=T),
  "mean" = apply(Full.divbyNull.Prot, 2, mean, na.rm=T),
  "fano" = apply(Full.divbyNull.Prot, 2, fano, na.rm=T)
))

# Normality tests --------------------------------------------------------------------------------------------
set.seed(21768134)
ccc = sample(richColors(ncol(Protz.NullNorm)+n)[-(1:n)] )

NormalityCheck = F
if (NormalityCheck) {
  pdfA4plot_on(pname = "Histograms.Normdata", rows = 4, cols=3)
  for (i in 1:ncol(Protz.NullNorm) ) {
    P =colnames(Protz.NullNorm)[i]
    g = na.omit.strip(Protz.NullNorm[,P])
    hist(g, xlab = P, breaks = 20, main = P, col=ccc[i], sub="Matched (mean, var) normal distribution overlaid")
    curve(dnorm(x, mean=mean(g), sd=sd(g)), lwd=2, lty=2,add=TRUE, yaxt="n")
  } #for
  pdfA4plot_off()
  
  pdfA4plot_on(pname = "QQ.plots.normality", rows = 4, cols=3)
  for (i in 1:ncol(Protz.NullNorm) ) {
    P =colnames(Protz.NullNorm)[i]
    qqPlot(Protz.NullNorm[,i], distribution="norm", lwd=.5, pch = 18, main = P
           , ylab = p0("Quantiles of ", P), xlab = "Quantiles of Normal Distribution"  )
  } #for
  pdfA4plot_off()
}


# Levene's test --------------------------------------------------------------------------------------------
DAT = get.oddoreven(Protz.NullNorm, odd = F)

ORD = colnames(DAT)
ORD2 = substr(ORD, 5, 15)
CN =colnames(DAT) = ORD2

PariwiseCombinations = data.table::CJ(CN,CN); idim(PariwiseCombinations)
pVal.l = Fstat.l = pVal = Fstat= matrix.fromNames(rowname_vec = ORD2, colname_vec = ORD2)

for (i in 1:NROW(PariwiseCombinations) ) {
  CPL = as.character(PariwiseCombinations[i,])
  #Melt data
  dataset <- melt(DAT[ ,CPL], na.rm = T)[2:3]
  if (CPL[1]==CPL[2]) { dataset$Var2 = as.factor(c(rep(1, NROW(dataset)/2), rep(2, NROW(dataset)/2) ))  } #if

# Levene's test
  #Compute test
  fit.l <- leveneTest(value ~ Var2, dataset )
  pVal.l[CPL[1], CPL[2]]   = fit.l$`Pr(>F)`[1]
  Fstat.l[CPL[1], CPL[2]]  = fit.l$`F value`[1]
}

annot_col.create.pheatmap.df(data = pVal, annot_df_per_column = Stat.prots[ORD2,c("var","mean")])
annot_col$var = c("white", "red")
annot_col$mean = c("white", "orange2")

mm =p0("Fig.S3D.Levene's test pairwise -log10(p-values)")
pheatmap::pheatmap(-log10(pVal.l), display_numbers = T, main = mm, annotation_col = annot, annotation_colors = annot_col, cutree_cols = 2)
wplot_save_this(mm, h=6)




AdditionalPlots= F
if (AdditionalPlots) {
  mm =p0("Levene's test pairwise p-values")
  pheatmap::pheatmap((pVal.l), display_numbers = T, number_format = "%.3f", main = mm, annotation_col = annot, annotation_colors = annot_col)
  wplot_save_this(mm, h=6)
  
  mm =p0("Levene's test pairwise F-statistics")
  pheatmap::pheatmap((Fstat.l), display_numbers = T, main = mm, annotation_col = annot, annotation_colors = annot_col)
  wplot_save_this(mm, h=6)
}

