######################################################################################################
# Kinetochore Regulation
######################################################################################################
# rm(list=ls(all.names = TRUE));  try(dev.off(), silent = T)
# source("~/GitHub/Kinetochore/Scripts/11.KinetochoreRegulation.R")

# Functions ------------------------
library(segmented)
library(gridExtra)

# Setup ------------------------

setup_MarkdownReports(OutDir = paste0(OutDirOrig,"11.KinetochoreRegulation"), scriptname = "11.KinetochoreRegulation.R")


# Data ------------------------
KinetochoreRegulation = as.data.frame(Final.Data.DF )

cnn = colnames(KinetochoreRegulation)
idx= seq(from = 1,to = ncol(KinetochoreRegulation), by=  2 )
protz = cnn[idx]

# Piecewise.Linear.regression ------------------------------------------------------------------------------------
Piecewise.Linear.regression = T
if (Piecewise.Linear.regression) {
  try.dev.off()
  Data.sec1 = Data.sec2 =list.fromNames(protz)
  Y.intercepts = BreakPoints = BreakPoints.y = Slopes =vec.fromNames(protz)

  pdfA4plot_on(pname = "Fig.2A.KinetochoreRegulation.piecewise", rows = 4, cols = 3)
  for (i in 1:l(idx)) {
    pr= protz[i]
    ix = idx[i]
    x = KinetochoreRegulation[, (ix+1):ix]
    cn=colnames(x)

    # linear modeling
    xm  =x; colnames(xm) = c("x","y")
    lin.mod <- lm( y~x, data = xm)
    segmented.mod <- segmented(obj = lin.mod, seg.Z = ~x)
    
    colnames(x)=c("k-fiber intensity/CENP-C", "Protein intensity/CENP-C")
    plot(x, pch=20, col="dodgerblue3", main = pr)
    points.segmented(segmented.mod)
    plot(segmented.mod, add = T)

    summary.segmented(segmented.mod)
    BreakPoints[pr] = iround(print.segmented(segmented.mod)$psi[2]) # Estimated Break-Point (x-dimension)
    Slopes[pr] = iround(print.segmented(segmented.mod)$coefficients["x"]) # Estimated coefficient (slope, for section 1)
    Y.intercepts[pr] = iround(segmented.mod$coefficients["(Intercept)"]) # Estimated Y axis intercept
    Data.sec1[[pr]] = x[,2][x[,1]<=BreakPoints[pr]]
    Data.sec2[[pr]] = x[,2][x[,1]>BreakPoints[pr]]
  } # for
  pdfA4plot_off()
} # if


VarianceInPhase2 = iround(unlapply(Data.sec2, var, na.rm = T))
MedianPhase2 = iround(unlapply(Data.sec2, median, na.rm = T))


# FULL data stat ---------------------------------------------
FULLexp = CENPC_norm_ls$FULL
idx_PRT = grep("BATCH|TUB|NN",colnames(FULLexp), invert = T, value = T)[1:7]
Medianz_Full.Attachment =  colMedians(FULLexp[,idx_PRT])


# NULL data stat ---------------------------------------------
NULLexp = CENPC_norm_ls$`NULL`
Medianz_NULL.empty =  colMedians(NULLexp[,idx_PRT])
Null2Full.median.decrease= Medianz_NULL.empty / Medianz_Full.Attachment

pname = "Fig.2C.Features.clustering.Zscore.with.Full"
Features = t(scale(cbind(
  'Null2Full.decrease [N/F]' = 	  Null2Full.median.decrease[str_split_fixed(protz, pattern = '\\.', n = 2)[,1]], 
  'BreakPoints [EXP]' = 				  BreakPoints[protz], 
  'Slopes [EXP]' = 					      Slopes[protz], 
  'VarianceInPhase2 [EXP]' = 		  VarianceInPhase2[protz], 
  'MedianPhase2 [EXP]' = 					MedianPhase2[protz], 
  'Medianz [Full]' = 					    Medianz_Full.Attachment[str_split_fixed(protz, pattern = '\\.', n = 2)[,1]]) 
) ) #

pheatmap::pheatmap(Features, treeheight_row = 20, treeheight_col = 20, main = pname, 
                   cluster_cols = T, cutree_cols = 2, clustering_method = "ward.D2")
wplot_save_this(plotname = pname, w = 6, h = 5)

# ----------------------------------------------------------------------------------------------------
