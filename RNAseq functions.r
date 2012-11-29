###
# RNA-seq functions
#source("C:/Users/jeremy/Documents/Code/useful_R.r")

write.delim = function(dat,File,sep='\t',row.names=FALSE,quote=FALSE,col.names=TRUE){
  write.table(dat,File,sep=sep,row.names=row.names,quote=quote,col.names=col.names)
}

### Generate all diff data
generateDiffData = function (cuff){
  diff.genes <<- diffData(genes(cuff))
  diff.isoforms <<- diffData(isoforms(cuff))
  diff.CDS <<- diffData(CDS(cuff))
  diff.TSS <<- diffData(TSS(cuff))
  F.genes <<- features(genes(cuff))
  F.isoforms <<- features(isoforms(cuff))
  F.CDS <<- features(CDS(cuff))
  F.TSS <<- features(TSS(cuff))
}



write.files = function(myData,fileBase){
  write.delim(cbind(F.genes,myData$genes),paste(fileBase,'_genes.txt',sep=''))
  write.delim(cbind(F.TSS,myData$TSS),paste(fileBase,'_TSS.txt',sep=''))
  write.delim(cbind(F.isoforms,myData$isoforms),paste(fileBase,'_isoforms.txt',sep=''))  
  write.delim(cbind(F.CDS,myData$CDS),paste(fileBase,'_CDS.txt',sep=''))
}  




### Function to get data between two samples
getPairData = function(myData, sample1, sample2){
  dataCol1 = cuffSampleNames[dataNames==sample1]
  dataCol2 = cuffSampleNames[dataNames==sample2]  
  if (length(dataCol1) < 1 | length(dataCol2) <1) stop(sprintf("Sample names must conform to %s",paste(dataNames,collapse=", ")))  
  
  # Retrieve data
  d1 = myData[(myData$sample_2 == dataCol1 & myData$sample_1 == dataCol2 ) | (myData$sample_1 == dataCol1 & myData$sample_2 == dataCol2),]
  # Rename the columns appropriately
  dimnames(d1)[[2]][dimnames(d1)[[2]]=='value_1'] <- dataNames[cuffSampleNames == d1$sample_1[1]]     
  dimnames(d1)[[2]][dimnames(d1)[[2]]=='value_2'] <- dataNames[cuffSampleNames == d1$sample_2[1]]
  # Return the data withouth the sample_1, sample_2
  d1[,!(dimnames(d1)[[2]] %in% c('sample_1','sample_2'))]
}

### Function to get the ALL the data for a comparison between 2 different samples
getAllPairData = function(sample1,sample2){
    list(genes = getPairData(diff.genes,sample1,sample2),isoforms = getPairData(diff.isoforms,sample1,sample2), CDS = getPairData(diff.CDS,sample1,sample2), TSS = getPairData(diff.TSS,sample1,sample2))
}

##################################################################
#### Functions to combine multiple lists
getSubset=function(x,SUBSET)
  subset(x,eval(parse(text=SUBSET)))
  
getMultipleSubsets = function(myData, SUBSET){
   
   ### First get genes
   myGenes = getSubset(myData$genes, SUBSET)
   geneIDs = myGenes$gene_id
   
   ### Now get TSS that don't appear in Genes
   allTSS = data.frame(gene_id = F.TSS$gene_id,myData$TSS)
   myTSS  = getSubset(allTSS[! (allTSS$gene_id %in% geneIDs),], SUBSET)
   myTSS2  = allTSS[allTSS$gene_id %in% myTSS$gene_id ,]
   
   # Get back all the gene Data   
   tssGenes = myData$genes[myData$genes$gene_id %in% unlist(strsplit(as.character(myTSS$gene_id),",")),]
   geneIDs = c(geneIDs,tssGenes$gene_id)
   
   ### Now get isoforms that don't appear in Genes or TSS
   allIsoforms = data.frame(gene_id = F.isoforms$gene_id,myData$isoforms)
   myIsoforms  = getSubset(allIsoforms[!(allIsoforms$gene_id %in% geneIDs),], SUBSET)
   myIsoforms2  = allIsoforms[allIsoforms$gene_id %in% myIsoforms$gene_id ,]
   
   # Now get back the next gene data
   isoGenes = myData$genes[myData$genes$gene_id %in% unlist(strsplit(as.character(myIsoforms$gene_id),",")),]
   geneIDs = c(geneIDs,isoGenes$gene_id)
  
   ### Now get CDS that don't appear elsewhere
   allCDS = data.frame(gene_id = F.CDS$gene_id,myData$CDS)
   myCDS  = getSubset(allCDS[! (allCDS$gene_id %in% geneIDs),], SUBSET)
   myCDS2  = allCDS[allCDS$gene_id %in% myCDS$gene_id ,]
   
   # Get back all the gene Data   
   cdsGenes = myData$genes[myData$genes$gene_id %in% unlist(strsplit(as.character(myCDS$gene_id),",")),]
   geneIDs = c(geneIDs,cdsGenes$gene_id)
     
     
   ## Combine all genes
   myGenes2 = rbind( cbind(Sig="1_gene",myGenes), cbind(Sig="2_TSS",tssGenes),cbind(Sig="3_isoform",isoGenes),cbind(Sig='4_CDS',cdsGenes)) 
   
   ## Merge with isoforms, tss, cds
   names(myTSS2) = paste("TSS",names(myTSS2),sep='_')
   names(myIsoforms2) = paste("isoform",names(myIsoforms2),sep='_')
   names(myCDS2) = paste("CDS",names(myCDS2),sep='_')
   
   allSig = merge(merge(merge(myGenes2,myTSS2,by.x = 2, by.y = 1,all=T), myIsoforms2,by=1,all=T),myCDS2,by=1,all=T)
   allSig = as.matrix(allSig[order(allSig$Sig,allSig$gene_id),])
   allSig[is.na(allSig)] <- ''
   allSig   
}  
  
  
### Go from fpkm data to a matrix ...
get_fpkm_matrix = function(genes,sampleNames=NULL){
  if(is.null(sampleNames)){
    sampleNames <- samples(genes) 
  }else{
    if(length(samples(genes)) != length(sampleNames))
      stop("Length of sample names must match the number of samples in the CuffData object")
  }
  
  f =fpkm(genes)
  myData = t(matrix(f$fpkm,nrow=length(sampleNames)))
  rownames(myData) = f$gene_id[!duplicated(f$gene_id)]
  dimnames(myData)[[2]] = sampleNames
  myData
}

log_or_NA = function(dat){
  dat2 = log10(dat)
  dat2[is.infinite(dat2)] = NA
  dat2
}

log_or_MIN = function(dat,MIN=NULL){   # You can specify a minimum value for 0 values to get after the log transform.  Defaults to the minimum non-0 value
  dat2 = log10(dat)
  if (is.null(MIN))
    MIN = min(dat2[!is.infinite(dat2)])
  #dat2[is.infinite(dat2)] = MIN
  dat2[dat2 < MIN] = MIN
  dat2
}

log2_or_MIN = function(dat,MIN=NULL){   # You can specify a minimum value for 0 values to get after the log transform.  Defaults to the minimum non-0 value
  dat2 = log2(dat)
  if (is.null(MIN))
    MIN = min(dat2[!is.infinite(dat2)])
  #dat2[is.infinite(dat2)] = MIN
  dat2[dat2 < MIN] = MIN
  dat2
}

heatmap1 = function(dat,main='',maxVal=NULL,Rowv=FALSE,margins=c(10,5),...){#samples=dimnames(dat)[[2]],geneNames=rownames(dat)){#
  require(gplots)

  if(is.null(maxVal)){
    RANGE = max(abs(range(dat,na.rm=T))) ;#*1.01, breaks=seq(-color_range,color_range))        
    breaks = seq(-RANGE,RANGE,length.out=20)
  }else
    breaks = seq(-maxVal,maxVal,length.out=20)
    h1=heatmap.2(dat,col=greenred,breaks=breaks,trace='n',Rowv=Rowv,margins=margins,Colv=FALSE,main=main,...)
}


norm_by_0 = function(dat,Hr0 = 1){
  zeroHr = dat[,Hr0]
  for(i in 1:ncol(dat))
    dat[,i] = dat[,i] - zeroHr
  dat
}

plot2Hists = function(X,main,xlab,xlim=c(0,50),breaks=seq(-0.5,1000.5,by=1),freqMax=seq(0,0.02,0.14),cols=c(rgb(0,0,1,alpha=0.4),rgb(1,0,0,alpha=0.4)),leg=NULL,...){
  par(new=F)
  mains = c(main,'');
  xaxt=c('n','s');ylab=c('Frequency','');xlab=c(xlab,'');
  ylim=c(0,max(freqMax))
  #print(ylim)
  cols2 = cols
  print(t.test(X[[1]],X[[2]]))
  a=lapply(X,function(x,main){
    a=hist(x,freq=F,breaks=breaks,xlim=xlim,ylim=ylim,col=cols[1],border=cols[1],main=mains[1],xaxt=xaxt[1],yaxt='n',ylab=ylab[1],xlab=xlab[1],cex.lab=1.5,cex.axis=1.2,cex.main=1.5,...);
    axis(2,at=freqMax,lab=freqMax*diff(breaks)[1],cex.axis=1.2)
    abline(v=median(x),col=cols[1],lwd=2);
    cols<<-cols[2];
    mains<<-mains[2];
    xaxt<<-xaxt[2];
    ylab<<-ylab[2];
    xlab<<-xlab[2];
    par(new=T)
    a
  }) 
  if(!is.null(leg)) legend('topright',fill=cols2,leg=leg)
  par(new=F)
  a
}





