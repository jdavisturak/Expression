## RNA-seq basic gene list analysis

#  1) Folder where you files live
setwd("C:/Users/Jeremy/Documents/Lab/Bryce/B_IgM_analysis/") 
# 2)  Location of helper file
source("C:/Users/jeremy/Documents/Code/RNAseq functions.r")
library(cummeRbund)

# 3) Read in data (slow the first time)
cuff=readCufflinks()
generateDiffData(cuff)  

# 4) Set sample names 
samples(cuff)$sample_name
## Option 1: names are currently 'q1', 'q2'....
dataNames = c("WT unstim", "WT 2h","WT 8h","WT 24h","Eko unstim","Eko 2h","Eko 8h","Eko 24h")
cuffSampleNames = paste('q',1:length(dataNames),sep='')

## Option 2:  it's already good..
#dataNames  = cuffSampleNames = samples(cuff)$sample_name

##################################################################################################

# 5) Pull out comparison of interest
DE_24 = getAllPairData('Eko 24h','WT unstim')
write.files(DE_24,'IgM_24')       # Write ALL files for this comparison (Genes, TSS, isoform, CDS)

# 6) Create a filter
DE_24_set1 = getMultipleSubsets(DE_24,"status=='OK' & q_value < 0.1")
write.delim(DE_24_set1,'IgM_24_set1.txt')  # write a file of combined data ('significant' genes only)

# 7) Heatmaps
# Get data all in one matrix
datNorm1 = norm_by_0(log2_or_MIN(get_fpkm_matrix(genes(cuff),c("WT unstim", "WT 2h","WT 8h","WT 24h","Eko unstim","Eko 2h","Eko 8h","Eko 24h"))))
#heatmap1(datNorm1[1:200,])

# Get subset of genes of interest
myGeneIDs = read.table("<<gene names list here >>")
heatmap1(datNorm1[toupper(rownames(datNorm1)) %in% toupper(myGeneIDs),])






