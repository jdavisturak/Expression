######################################################################################
######################################################################################
### Simple_contrasts.r  
##  Jeremy Davis-Turak  2/22/2012
##  This script contains functions that help in setting up microarray analysis of 
##  expression data from various timepoints and stimuli.  
##  
######################################################################################
######################################################################################

## Normalization functions
saveWithGeneNames = function(data,genes,file)
	write.table(cbind(genes,data),file,sep="\t",quote=F,row.names=F)
subtract0 = function(row) row-row[1]

### Function to do Venn Diagram-type data
overlap2or3 = function(outputName, A, B, C=NULL,N=c('A','B','C'), column='SYMBOL' ){
	
	# Get rid of duplicate stuff, keep only names/probes
	A_nam = as.character(A[,column])
	A_names = A_nam[!duplicated(A_nam)]
	B_nam = as.character(B[,column])
	B_names = B_nam[!duplicated(B_nam)]

	Two_AB = intersect(A_names,B_names)
	Two_B = setdiff(B_names,A_names)
	Two_A = setdiff(A_names,B_names)
	
	if(!is.null(C)){
	
		C_nam = as.character(C[,column])
		C_names = C_nam[!duplicated(C_nam)]
	
		Three_A = setdiff(Two_A,C_names)
		Three_B = setdiff(Two_B,C_names)
		
		Three_AC = intersect(Two_A,C_names)
		Three_BC = intersect(Two_B,C_names)
		
		Three_ABC = intersect(Two_AB,C_names)
		Three_AB = setdiff(Two_AB,C_names)
		
		Three_C = setdiff(C_names, c(Two_AB,Two_B,Two_A))
		
		## Rename, for laziness
		Two_A = Three_A
		Two_B = Three_B
		Two_AB = Three_AB
		
	}
	
	### Make summary table:
	# A, B, AB, C, AC, BC, ABC
	categories = c(N[1],N[2],paste(N[1],N[2],sep= " AND "))
	nums = c(length(Two_A),length(Two_B),length(Two_AB))	
	#browser()
	data = rbind(
		cbind(A[A_nam %in% Two_A,],categories=categories[1]), 
		cbind(B[B_nam %in% Two_B,],categories=categories[2]), 
		cbind(B[B_nam %in% Two_AB,],categories=categories[3]))
	
	if(!is.null(C)) {
		categories = c(categories,N[3],paste(N[1],N[3],sep= " AND "),paste(N[2],N[3],sep= " AND "),paste(N[1],N[2],N[3],sep= " AND "))
		nums = c(nums,length(Three_C),length(Three_AC),length(Three_BC),length(Three_ABC))		
		#browser()
		data = rbind(data,
			cbind(C,categories=categories[4])[C_nam %in% Three_C,],
			cbind(C,categories=categories[5])[C_nam %in% Three_AC,],
			cbind(C,categories=categories[6])[C_nam %in% Three_BC,],
			cbind(C,categories=categories[7])[C_nam %in% Three_ABC,])		
	}
	
	
	stats = cbind(categories,nums)
	print(stats)
	
	write.table(data,file=paste(outputName,'_list.txt',  sep=""),sep="\t",quote=F,row.names=F)
	write.table(stats,file=paste(outputName,'_stats.txt',  sep=""),sep="\t",quote=F,row.names=F,col.names=F)
	
	
	list(stats=stats,list=data)
}

### make design
makeDesign = function(sample_names){
	design <- model.matrix(~0 + factor(sample_names))
	colnames(design) = levels(factor(sample_names))
	rownames(design) = sample_names
}

### Function to add on to the name of a sample string (except if it ends with 0)
concat = function(string,add){
  sub("([^0])$",sprintf("\\1%s",add),string)
}

### Function to add on to the name of EVERY sample in a contrast string (except if it ends with 0)
# Also, each name must be at least 2 characters long
concatC = function(string,add){
  gsub("([A-Za-z0-9_][A-Za-z1-9_])([^A-Za-z0-9_]|$)",sprintf("\\1%s\\2",add),string)
}

#### Function to replace the names of a matrix
replaceNames = function(mat,levels,contrasts,add){
  dimnames(mat)[[1]] = concat(levels,add)
  dimnames(mat)[[2]] = concatC(contrasts,add)
  mat
}

### Function to start with a matrix, and add another one after the diagonal...
addMatrix = function(mat1,mat2){
  if(is.null(mat1)) return(mat2);
  #rbind( cbind(mat1,matrix(0,nr=nrow(mat1),nc=ncol(mat2))), cbind(matrix(0,nr=nrow(mat2),nc=ncol(mat1)), mat2))
  out = cbind(rbind(mat1, matrix(0,nr=nrow(mat2),nc=ncol(mat1))),  rbind(matrix(0,nr=nrow(mat1),nc=ncol(mat2)), mat2))
  rownames(out)= c(rownames(mat1),rownames(mat2))
  out
}

### Function to get all possible combinations of stimuli and times
makeNames = function(stims,times)
  paste("_",as.vector(t(outer(stims,times,paste,sep="_"))),sep="")
  
  
  
### Function to combine the rows of a matrix (get rid of duplicate row names)
consolidateRows = function(mat){
  # Find duplicated row names
  Names = rownames(mat)
  dup.names = unique(Names[duplicated(Names)])
  inDups = Names %in% dup.names

  # Keep non-dupiplicated rows
  mat1 = subset(mat,!inDups)
  
  # For duplicated rows, perform a sort of 'OR' calculation to see which contrasts it's in...
  for(nam in dup.names){
    theseRows = mat[Names==nam,]
    # For each contrast (column), get the unique numbers (should only be 0 and 1, or 0 and -1, if there is more than one unique value)
    outputVals = numeric(ncol(mat))
    for(i in 1:ncol(mat)){
      these.vals = unique(theseRows[,i])
      realVal =  switch(length(these.vals),these.vals, these.vals[these.vals!=0]) # if len=1,return val. if len=2, return the non-zero. otherwise, NULL.
      if(length(realVal) != 1)
        stop(sprintf("Mal-formed design matrix. Row %s has more than 1 non-zero value in contrast %s",nam,colnames(mat)[i]))
        
      outputVals[i] = realVal  
    }            
    mat1 = rbind(mat1,outputVals)   # Add this row back
    rownames(mat1)[nrow(mat1)] = nam  # Give it the correct row name
  }
  mat1
}


  
### OK, given a matrix, stim and time names, output a new matrix
makeMatrix = function(mat1,stims,times,design){
  mat2 = NULL
  suffixes = makeNames(stims,times);
  for(i in 1:length(suffixes)){
    temp.mat = replaceNames(mat1,rownames(mat1),colnames(mat1),suffixes[i])
    temp.mat = temp.mat[, !(colnames(temp.mat) %in% colnames(mat2)),drop=F] # get rid of duplicate contrast names
    mat2 = addMatrix(mat2,temp.mat)
  }
  mat2 = consolidateRows(mat2)  # Get rid of duplicate rows
  #mat2 = mat2[,apply()] # do I need to get rid of contrasts that won't work?
  cleanContrasts(mat2,design)
}


### Function to clean up contrast matrix
cleanContrasts = function(mat,design){
	
	## Make sure all coefficients of design are respresented.
	missing.designs = setdiff(colnames(design),rownames(mat))
	matBlank = matrix(0, nr=length(missing.designs), nc=ncol(mat))
	rownames(matBlank) = missing.designs		
	mat1 = rbind(mat,matBlank)
	
	## Which rows of the contrast matrix are in fact not in the original data?
	spurious.rows = setdiff(rownames(mat1),colnames(design))
	badContrasts = NULL
	
	for( s in spurious.rows){
		badContrasts = c(badContrasts, which(mat1[s,] != 0)) # Bad				
	}
	
	# Get rid of those columns
	if(length(badContrasts) > 0)
		mat1 = mat1[,-unique(badContrasts)]
	# Get rid of those rows
	mat1 = mat1[!(rownames(mat1) %in% spurious.rows),]
	
	## Get rid of contrasts that are all 0
	mat1 = mat1[,!apply(mat1, 2, function(x) all(x==0))]
	
	## Re-order columns of the contrast matrix
	
			
}


computeContrasts=function(log2data,design,contrast.matrix){
	design = design[,order(colnames(design))]
	contrast.matrix = contrast.matrix[order(rownames(contrast.matrix)),]
	contrasts.fit(lmFit(log2data, design), contrast.matrix)
}

saveWithGeneNames = function(data,genes,file)
	write.table(cbind(genes,data),file,sep="\t",quote=F,row.names=F)


saveAndReturnGenes = function(masterGenes,newIndex,name){
	out = masterGenes[newIndex, ]
	genes = normData[newIndex, 1:6]
	saveWithGeneNames(out, genes,name)
	out
}
#########################################################################
#### Sample code

## Read in targets
## Read in data

## Normalize data?
# (either by time 0 of WT, or time 0 of specific genotypes)

## set up a contrast matrix
# Example: 
#           WT-WT_0   p50_hyper_induced  p50_hyper_basal   
            
# WT            1           -1
# WT_0         -1            0
# p50KO         0            1
# p50KO_0       0            0
# IFNARKO       0            0
# IFNARKO_0     0            0
# DOUBLEKO      0            0
# DOUBLEKO_0    0            0


#makeMatrix(cm,c("LPS","CPG"),1:3)
