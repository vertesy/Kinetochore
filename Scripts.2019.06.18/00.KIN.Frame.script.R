######################################################################################################
# 00.KIN.Frame.script.R
######################################################################################################
rm(list=ls(all.names = TRUE));  try(dev.off(), silent = T)
# source("~/GitHub/Kinetochore/Scripts.2019.06.18/00.KIN.Frame.script.R")

# Functions ------------------------
# If compatibity issues arise, use MarkdownReports v2.9.5 at: https://github.com/vertesy/MarkdownReports/releases/tag/v2.9.5
# install.packages("devtools"); # If you don't have it
source ('/Users/abel.vertesy/GitHub/Kinetochore/Scripts.2019.06.18/CodeAndRoll.v2018.07.03.R')
require("ggplot2")
require("stringr")
require("devtools")
# devtools::install_github(repo = "vertesy/MarkdownReports.v2.9.5")
require(MarkdownReports)
# require("MarkdownReportsDev")


# Setup ------------------------
ScriptDir = "~/GitHub/Kinetochore/Scripts.2019.06.18/"
InputDirOrig = "~/GitHub/Kinetochore/Data/"
OutDirOrig = "~/Google_Drive/Avano/KIN/Analysis.v2019.06.19.v4/"
setup_MarkdownReports(OutDir = OutDirOrig, scriptname = "00.KIN.Frame.script.R")

# Parameters ------------------------

usenew=T # use new data d2017.08.23
SHEETZ = c( 'EXP', 'NULL', 'FULL')

# 01. Normalisation Steps on the raw data ------------------------
source("~/GitHub/Kinetochore/Scripts.2019.06.18/01.Raw.Data.Analysis.R")

# 05 Null normalised ------------------------
source("~/GitHub/Kinetochore/Scripts.2019.06.18/05.Fig3.Tests.Normality.Variance.R")

# Experiment Comparison  ------------------------
source("~/GitHub/Kinetochore/Scripts.2019.06.18/10.Experiment.Comparison.R")

# 11.KinetochoreRegulation.R ------------------------
source("~/GitHub/Kinetochore/Scripts.2019.06.18/11.KinetochoreRegulation.R")


#  ------------------------
