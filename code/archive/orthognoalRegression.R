#!/usr/bin/env Rscript
#======================================================================================#
#  This code produces a graph depicting the orthoganol regression
#    between two species. Orthogonal regression is an approach that minimizes
#    the perpendicular distances between the data points (unlike traditional 
#    linear regressions that only minimize the sum of squared vertical 
#    distances on the y-axis). 
#
#  Regressions are performed on geometric mean values for each chemical. This
#    code also plots the raw experimental values, allowing the reader to see the
#    full experimental variability for each chemical.
#
#  This code recreates Figure 3 in Connors et al. 2022. ETC. 41(1)134-147. 
#
#  Code requires a data.table with columns named "CAS", "Latin.name", and "Effect.value"
#======================================================================================#
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# source functions
source("code/functions.R")

# load data
d1.c <- data.table::fread("data/rodentNematode.csv")

#===========================================================#
# New code for orthogonal regression
#===========================================================#


#===========================================================#
# Original ploting code working from the unkown rat.csv file
#===========================================================#
# load data
# d1.c <- data.table::fread("data/rat.csv") # this file was not sent to me TAC. Tried to recreate
d1.c <- data.table::fread("data/rodentNematode.csv")

# pdf(file="Paper Figures Test.pdf",height=8,width=8)
doTitles<-FALSE
#pdf(file="Figures/MainFigures99.9.pdf",height=8,width=8,paper="USr")
par2start <- par()
par(pty='s',omi=c(0,0,0,0),mai=c(.8,.8,.1,.1))
SpeciesList<-c("Rodent","Caenorhabditis elegans")
ChemList<-unique(c(d1.c[Latin.name=="Rodent", CAS], d1.c[Latin.name=="Caenorhabditis elegans", CAS]))
isJava<-(substring(names(dev.cur()),1,4)=="java")
pageno <- 0

lineColor1 <- lineColor2 <- rgb(252, 141, 89,maxColorValue=256)
refColor <- "cyan"
pointsCOL <- "red"
grayCEX <- 0.8

dm <- data.table::as.data.table(d1.c[Latin.name=="Rodent", list(CAS, Effect.value)])
##  plotting algorithm expects data in units ug/L. My data is in mg/L, adjusting units manually.
dm[, Effect.value :=  Effect.value * 1000]
dm <- as.data.frame(dm)


cd <- data.table::as.data.table(d1.c[Latin.name=="Caenorhabditis elegans", list(CAS, Effect.value)])
##  plotting algorithm expects data in units ug/L. My data is in mg/L, adjusting units manually.
cd[, Effect.value :=Effect.value * 1000]
cd <- as.data.frame(cd)

plotOutput <- plotsForPaper(
  ySpecies=expression(paste(italic("Rodent"), " toxicity", phantom(0),LD50,phantom(0), "(mg/kg)")),
  xSpecies=expression(paste(italic("Caenorhabditis elegans"), " toxicity",phantom(0),EC10,phantom(0),"(",mu,"g/L)")),
  chemList = ChemList,
  yData = dm,
  yVar = "Effect.value",
  xData = cd,
  xVar = "Effect.value",
  mainTitle = "",
  pageText = paste("Em vs El, pg.",pageno),
  meansOnly=FALSE,
  refCOL="#80FFFFE6",
  regCOL=lineColor1,
  confidenceAlpha=5
)

# addRegInfo(plotOutput,corr=T,eqnCOL=lineColor1)
dim(plotOutput[[1]])
if(doTitles)mtext(side=3,line=-1,adj=0,text="D.magna and C.dubia Chronic Toxicity")
# dev.off()


