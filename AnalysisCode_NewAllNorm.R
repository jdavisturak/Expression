library(gplots)
library(limma)
setwd("/Users/Kristyn/Desktop/Microarray!!/Kristyn Microarray Analysis/")
source("Simple_contrasts.r")

normData = read.delim("hoffmann+IJLMO_norm_KEF data.txt",stringsAsFactors=F)
targets = read.delim("targets_120215.txt",stringsAsFactors=F)
targets$GenoType = gsub("/","_",targets$GenoType)

##### Log 2 transform
rawNormData = normData[,-(1:6)]     
log2NormData = log2(rawNormData)

## Here is the place to change targets

### Make a 'design' table where rows are individual samples, and the columns represent unique categories (i.e. group the replicates)
sample_names = paste(targets$GenoType,targets$Batch,targets$Stimulus,targets$Hour,sep="_")
sample_labels = paste(targets[,1],sample_names)

#design = makeDesign(sample_names)
design <- model.matrix(~0 + factor(sample_names))
	colnames(design) = levels(factor(sample_names))
	rownames(design) = sample_names
	
## Tabulate all sample types
apply(design,2,function(x)length(which(x==1)))->a
a[a>1]


# #### Sample checking  (outliers)            
# # Correlate every sample ...    
# IAC = cor(rawNormData)               
# heatmap.2(IAC,trace='n',margins=c(10,10),labCol= sample_labels ,labRow= sample_labels               )

# ## Get rid of number 23
# outliers = 23
# h1=heatmap.2(1-IAC[-outliers,-outliers],trace='n',margins=c(10,10),labCol= sample_labels[-outliers] ,labRow= sample_labels[-outliers])
# sample_names[h1$rowInd[1:14]]

# ### Call that one sample 'bad'
# sample_names[23] = paste(sample_names[23],"bad",sep="_")
# design = makeDesign(sample_names)


##############################################################################################################################
## THIS IS THE IMPORTANT PART TO EDIT

### Set up the contrast matrix
# Rows are different sample categories, columns are comparison 
temp_row_names = sort(c(unique(paste(targets$GenoType,targets$Batch,sep="_")),unique(sample_names[regexpr("0$",sample_names)>0])))
temp_row_names
 # [1] "IFNAR_IRF3dko"      "IFNAR_IRF3dko_NS_0" "IFNAR_p50dko"       "IFNAR_p50dko_NS_0"  "IFNARko"            "IFNARko_NS_0"      
 # [7] "p50ko"              "p50ko_NS_0"         "WT"                 "WT_NS_0"                   

write.table(temp_row_names,file='temp_row_namesCompiled.txt',sep="\t",quote=F)

##############################################################################################################################
# Can repeat this part of the code
contrast.matrix = as.matrix(read.delim('temp_row_names2.txt',stringsAsFactors=F)[,-1,drop=F])
rownames(contrast.matrix) = temp_row_names
contrast.matrix2=makeMatrix(contrast.matrix,c("LPS", "CPG", "IFN", "PIC", "TNF"),c(1,3,8,24),design)

## Get Data
fitData = computeContrasts(log2NormData,design,contrast.matrix2)
temp_logRatios = fitData$coefficients
colnames(temp_logRatios )
##############################################################################################################################

#logRatios = fitData$coefficients[,c(1,2,3)]
logRatios = fitData$coefficients
#write.table(cbind(normData[,1:6], logRatios),file='test.txt',sep="\t",quote=F)

########################################## Pull out genes KF&CC WT, LPS, 1,3,8
WT_LPS_columns_aa = c("WT_LPS_1.WT_NS_0", "WT_LPS_3.WT_NS_0", "WT_LPS_8.WT_NS_0")
WT_LPS_columns_cc = c("WT_cc_LPS_1.WT_cc_NS_0", "WT_cc_LPS_3.WT_cc_NS_0", "WT_cc_LPS_8.WT_cc_NS_0")

logRatioForInduced = log2(3)
WT_LPS_up = apply(logRatios[,WT_LPS_columns_aa],1,max) > logRatioForInduced
table(WT_LPS_up)

WT_LPS_up_genes  = saveAndReturnGenes(logRatios[,WT_LPS_columns_aa],WT_LPS_up,sprintf("KF_LPS_WT_induced_LogRatio%.2f.txt", logRatioForInduced))
WT_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(WT_LPS_columns_aa, WT_LPS_columns_cc)],WT_LPS_up,sprintf("KF_LPS_WT_induced_LogRatio%.2fAll.txt", logRatioForInduced))

### Pull out genes CC, WT, LPS, 1,3,8
logRatioForInduced_cc = log2(3)
WT_LPS_up_cc = apply(logRatios[,WT_LPS_columns_cc],1,max) > logRatioForInduced
table(WT_LPS_up_cc)

CC_WT_LPS_up_genes  = saveAndReturnGenes(logRatios[,WT_LPS_columns_cc],WT_LPS_up_cc,sprintf("CC_LPS_WT_induced_LogRatio%.2f.txt", logRatioForInduced))
CC_WT_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(WT_LPS_columns_aa, WT_LPS_columns_cc)],WT_LPS_up_cc,sprintf("CC_LPS_WT_induced_LogRatio%.2fAll.txt", logRatioForInduced))

######
# dat1 = read.delim('KristynData_LPS_WT_induced_LogRatio1.00All.txt',stringsAsFactors=F)
# dat2 = read.delim('CC_Data_LPS_WT_induced_LogRatio1.00All.txt',stringsAsFactors=F)
# a=overlap2or3('test1', dat1, dat2, N=c('lps','lps(cc)'))

dat1 = read.delim('KF_LPS_WT_induced_LogRatio1.58All.txt',stringsAsFactors=F)
dat2 = read.delim('CC_LPS_WT_induced_LogRatio1.58All.txt',stringsAsFactors=F)
a=overlap2or3('wt_lps_venn', dat1, dat2, N=c('lpsKF','lpsCC'))

##################################################################################################################
### Pull out genes KF&CC p50, LPS, 1,3,8, norm to WT0
p50_LPS_columns_aa = c(6,8,24,37)
p50_LPS_columns_cc = c(3,5,22,35)

logRatioForInduced = log2(3)
p50_LPS_up = apply(logRatios[,p50_LPS_columns_aa],1,max) > logRatioForInduced
table(p50_LPS_up)

p50_LPS_up_genes  = saveAndReturnGenes(logRatios[,p50_LPS_columns_aa],p50_LPS_up,sprintf("KF_LPS_p50_induced_LogRatio%.2f.txt", logRatioForInduced))
p50_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(p50_LPS_columns_aa, p50_LPS_columns_cc)],p50_LPS_up,sprintf("KF_LPS_p50_induced_LogRatio%.2fAll.txt", logRatioForInduced))

### Pull out genes CC, p50, LPS, 1,3,8
logRatioForInduced_cc = log2(3)
p50_LPS_up_cc = apply(logRatios[,p50_LPS_columns_cc],1,max) > logRatioForInduced
table(p50_LPS_up_cc)

CC_p50_LPS_up_genes  = saveAndReturnGenes(logRatios[,p50_LPS_columns_cc],p50_LPS_up_cc,sprintf("CC_LPS_p50_induced_LogRatio%.2f.txt", logRatioForInduced))
CC_p50_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(p50_LPS_columns_aa, p50_LPS_columns_cc)],p50_LPS_up_cc,sprintf("CC_LPS_p50_induced_LogRatio%.2fAll.txt", logRatioForInduced))

dat1 = read.delim('KF_LPS_p50_induced_LogRatio1.58All.txt',stringsAsFactors=F)
dat2 = read.delim('CC_LPS_p50_induced_LogRatio1.58All.txt',stringsAsFactors=F)
a=overlap2or3('p50_lps_venn', dat1, dat2, N=c('lpsKF(p50)','lpsCC(p50)'))

##################################################################################################################
### Pull out genes KF&CC IFNAR, LPS, 1,3,8, norm to IFNAR0
IFNAR_LPS_columns_aa = c(16,30,43)
IFNAR_LPS_columns_cc = c(17,31,44)

logRatioForInduced = log2(3)
IFNAR_LPS_up = apply(logRatios[,IFNAR_LPS_columns_aa],1,max) > logRatioForInduced
table(IFNAR_LPS_up)

IFNAR_LPS_up_genes  = saveAndReturnGenes(logRatios[,IFNAR_LPS_columns_aa],IFNAR_LPS_up,sprintf("KF_LPS_IFNAR_induced_LogRatio%.2f.txt", logRatioForInduced))
IFNAR_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(IFNAR_LPS_columns_aa, IFNAR_LPS_columns_cc)],IFNAR_LPS_up,sprintf("KF_LPS_IFNAR_induced_LogRatio%.2fAll.txt", logRatioForInduced))

### Pull out genes CC, IFNAR, LPS, 1,3,8
logRatioForInduced_cc = log2(3)
IFNAR_LPS_up_cc = apply(logRatios[, IFNAR_LPS_columns_cc],1,max) > logRatioForInduced
table(IFNAR_LPS_up_cc)

CC_IFNAR_LPS_up_genes  = saveAndReturnGenes(logRatios[, IFNAR_LPS_columns_cc],IFNAR_LPS_up_cc,sprintf("CC_LPS_IFNAR_induced_LogRatio%.2f.txt", logRatioForInduced))
CC_IFNAR_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(IFNAR_LPS_columns_aa, IFNAR_LPS_columns_cc)],IFNAR_LPS_up_cc,sprintf("CC_LPS_IFNAR_induced_LogRatio%.2fAll.txt", logRatioForInduced))

dat1 = read.delim('KF_LPS_IFNAR_induced_LogRatio1.58All.txt',stringsAsFactors=F)
dat2 = read.delim('CC_LPS_IFNAR_induced_LogRatio1.58All.txt',stringsAsFactors=F)
a=overlap2or3('IFNAR_lps_venn', dat1, dat2, N=c('lpsKF(ifnar)','lpsCC(ifnar)'))


### Pull out HYPEREXPREssed genes KF&CC p50, LPS, 1,3,8, 
p50_LPS_columns_aa = c(2,9,14)
p50_LPS_columns_cc = c(1,8,13)

logRatioForInduced = log2(2)
p50_LPS_up = apply(logRatios[,p50_LPS_columns_aa],1,max) > logRatioForInduced
table(p50_LPS_up)

p50_LPS_up_genes  = saveAndReturnGenes(logRatios[,p50_LPS_columns_aa],p50_LPS_up,sprintf("KF_LPS_p50_hyper_LogRatio%.2f.txt", logRatioForInduced))
p50_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(p50_LPS_columns_aa, p50_LPS_columns_cc)],p50_LPS_up,sprintf("KF_LPS_p50_hyper_LogRatio%.2fAll.txt", logRatioForInduced))

### Pull out genes CC, p50, LPS, 1,3,8
logRatioForInduced_cc = log2(3)
p50_LPS_up_cc = apply(logRatios[,p50_LPS_columns_cc],1,max) > logRatioForInduced
table(p50_LPS_up_cc)

CC_p50_LPS_up_genes  = saveAndReturnGenes(logRatios[,p50_LPS_columns_cc],p50_LPS_up_cc,sprintf("CC_LPS_p50_hyper_LogRatio%.2f.txt", logRatioForInduced))
CC_p50_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(p50_LPS_columns_aa, p50_LPS_columns_cc)],p50_LPS_up_cc,sprintf("CC_LPS_p50_hyper_LogRatio%.2fAll.txt", logRatioForInduced))

dat1 = read.delim('KF_LPS_p50_hyper_LogRatio1.58All.txt',stringsAsFactors=F)
dat2 = read.delim('CC_LPS_p50_hyper_LogRatio1.58All.txt',stringsAsFactors=F)
a=overlap2or3('p50hyper_lps_venn', dat1, dat2, N=c('lpsKF(p50hyper)','lpsCC(p50hyper)'))


##################################################################################################################
### Pull out HYPOexpressed genes IFNAR/IRFS, LPS, 1,3,8, from IFNAR
IRF3_LPS_columns_aa = c(5,12,17)
Ip50_LPS_columns_aa = c(3,10,15)

logRatioForInduced = log2(3)
IRF3_LPS_up = apply(logRatios[,IRF3_LPS_columns_aa],1,max) > logRatioForInduced
table(IRF3_LPS_up)

IRF3_LPS_up_genes  = saveAndReturnGenes(logRatios[,IRF3_LPS_columns_aa],IRF3_LPS_up,sprintf("KF_LPS_IRF3_hypo_LogRatio%.2f.txt", logRatioForInduced))
IRF3_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(IRF3_LPS_columns_aa, Ip50_LPS_columns_aa)],IRF3_LPS_up,sprintf("KF_LPS_IRF3_hypo_LogRatio%.2fAll.txt", logRatioForInduced))

logRatioForInduced = log2(2.3)
Ip50_LPS_up = apply(logRatios[,Ip50_LPS_columns_aa],1,max) > logRatioForInduced
table(Ip50_LPS_up)

Ip50_LPS_up_genes  = saveAndReturnGenes(logRatios[,Ip50_LPS_columns_aa],Ip50_LPS_up,sprintf("KF_LPS_Ip50_hyper_LogRatio%.2f.txt", logRatioForInduced))
Ip50_LPS_up_genesAll  = saveAndReturnGenes(logRatios[,c(IRF3_LPS_columns_aa, Ip50_LPS_columns_aa)],Ip50_LPS_up,sprintf("KF_LPS_Ip50_hyper_LogRatio%.2fAll.txt", logRatioForInduced))

##Venn for IRF3 and IFNAR/p50
dat1 = read.delim('KF_LPS_IRF3_hypo_LogRatio1.58All.txt',stringsAsFactors=F)
dat2 = read.delim('KF_LPS_Ip50_hyper_LogRatio1.00All.txt',stringsAsFactors=F)
a=overlap2or3('IFNARp50_IRF3_lps_venn2', dat1, dat2, N=c('IRF3','IFNAR/p50'))














                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                