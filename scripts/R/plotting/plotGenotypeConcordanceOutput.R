#title           :processGenotypeConcordanceOutput.R
#description     :This script will analysis and plot the Output of GATK GenotypeConcordance
#author		 :Maximilian Hastreiter
#date            :20141210
#version         :0.1
#usage		 :Rscript processGenotypeConcordanceOutput.R <INFILE>
#==============================================================================

library(VennDiagram)
library(beeswarm)

#Get Infile
args <- commandArgs(TRUE)
infile = args[1]
eval   = args[2]
comp   = args[3]

pdf(paste(infile,"_",comp,"_SummaryPlots.pdf",sep=""))

SLS = read.table(paste(infile,"_SiteLevelSummaryStatistics",sep=""),head=T)
draw.pairwise.venn(SLS[5]+SLS[1],SLS[6]+SLS[1],as.numeric(SLS[1]),c(eval, comp),scaled=T)


#PGC = read.table(paste(infile,"_ProportionsGenotypesComp",sep=""),head=T, row.names=1)
#CC = read.table(paste(infile,"_ComparisonCounts",sep=""),head=T, row.names=1)
#PGE = read.table(paste(infile,"_ProportionsGenotypesEval",sep=""),head=T, row.names=1)

PSS = read.table(paste(infile,"_PerSampleSummaryStatistics",sep=""),head=T, row.names=1)
beeswarm(PSS, labels=c("NRS","NRD","OGC"),main="PerSampleSummaryStatistics")
boxplot(PSS,add=T,names=c("","",""))
dev.off()
