# This module contains functions to do the following:
# * Load TPM or Count data from batches and merge into a single giant sparse data matrix
# * Load QC data from batches and merge into a single meta-data file
# * convert from transcript-level counts or TPM to gene-level counts or TPM
# Typically one should run data_assemble.py first
# Primary function to use is load_all_samples
# Input required is simply the list of batch identifiers, such as: c("tumor1","tumor2")
 
library(modules)
import_package("Matrix",attach=TRUE)
parallel<-import_package("parallel")
sqldf<-import_package("sqldf")
plyr<-import_package("plyr")
biomaRt<-import_package("biomaRt") #converts transcript ID to gene ID
#fns<-import("../util/functions")

load_sample<-function(sample_name,transcript2gene_map=NULL,new_cell_id=TRUE,tpm=TRUE,folder="data"){
  #set tpm=FALSE to get counts
  #if transcript2gene_map is provided, gene-level values returned, 
  #otherwise, returns transcript-level
  # new_cell_id=TRUE means prepend sample name before cell name to make unique ID
  # new_cell_id=FALSE means use existing cell name as unique ID
  #transcript2gene_map should be a matrix with two columns:
  #first column contains transcript IDs
  #second column contains corresponding gene IDs
  f<-file.path(folder,sample_name,"abundance_combined.csv")
  d<-read.csv(f)
  d$sample<-factor(rep(sample_name,nrow(d)))
  d$length<-NULL #save memory space
  d$eff_length<-NULL
  if(new_cell_id==TRUE){
    d$cell_id<-factor(paste0(sample_name,"-",d$cell))
  } else {
    d$cell_id<-d$cell
  }
  if(tpm==FALSE){
    d$tpm<-NULL
  } else {
    d$est_counts<-NULL
  }
  id_map<-transcript2gene_map
  if(is.null(id_map)){
    colnames(d)[colnames(d)=="target_id"]<-"ensembl_transcript_id"
    return(d)
  }
  #this part run only if id_map is provided (ie, want to convert transcript ID to gene ID)
  colnames(id_map)<-c("transcript_id","gene_id")
  colnames(d)[colnames(d)=="target_id"]<-"transcript_id"
  d<-plyr$join(d,id_map,by="transcript_id")
  #some transcripts may not map to genes, so use transcript ID for these?
  #sqldf("select count(*) from d where gene_id is NULL") #typically 0.2% of transcripts
  #nagene<-is.na(d$gene_id)
  #d$gene_id[nagene]<-d$transcript_id[nagene]
  d<-d[complete.cases(d),] #currently unmappable transcripts are discarded!
  d$transcript_id<-NULL
  #note ddply is too slow, so use tapply or sqldf instead
  if(tpm==FALSE){
    #d<-ddply(d,c("cell_id","gene_id"),summarise,est_counts=sum(est_counts))
    d<-sqldf$sqldf("select cell_id,sample,gene_id,sum(est_counts) as est_counts from d group by cell_id,gene_id")
  } else {
    d<-sqldf$sqldf("select cell_id,sample,gene_id,sum(tpm) as tpm from d group by cell_id,gene_id")
  }
  return(d)
}

coo2matrix<-function(row_ids,col_ids,vals){
  # provide a list of row ids, a list of column ids, and a list of values
  # creates sparse Matrix
  rows<-factor(row_ids)
  cols<-factor(col_ids)
  mat <- Matrix(0.0, nrow=nlevels(rows), ncol=nlevels(cols))
  rownames(mat)<-levels(rows)
  colnames(mat)<-levels(cols)
  # Notice that you need to specify zero values to make it sparse.
  mat[cbind(rows,cols)] <- vals
  return(mat)
}

qc2meta<-function(sample_name,new_cell_id=TRUE,folder="data"){
  #read in a qc file, output a modified data frame with the batch label added
  f<-file.path(folder,sample_name,"kallisto_qc.csv")
  qcdat<-read.csv(f)
  colnames(qcdat)[1]<-"cell" #assume first column is the cell id
  qcdat$sample<-factor(rep(sample_name,nrow(qcdat)))
  if(new_cell_id==TRUE){
    qcdat$cell_id<-factor(paste0(sample_name,"-",qcdat$cell)) #useful for merging
  } else {
    qcdat$cell_id<-factor(qcdat$cell)
  }
  return(qcdat)
}

load_all_samples<-function(sample_list,new_cell_id=TRUE,transcripts2genes=TRUE,tpm=TRUE,biomart_dataset="hsapiens_gene_ensembl",biomart_cache=TRUE){
  # provide list of batch identifiers "sample_list" 
  # transcripts2genes: flag indicating whether to convert to gene counts from transcript counts input
  # new_cell_id=TRUE means prepend sample name before cell name to make unique ID
  # new_cell_id=FALSE means use existing cell name as unique ID
  # biomart_database="hsapiens_gene_ensembl" for human, "mmusculus_gene_ensembl" for mouse
  # biomart_cache forces using June 2015 database to avoid server problems with recent updates
  # returns a list with two components
  # "m" is a sparse matrix whose rows are genes/transcripts and cols are cells
  # "meta" is metadata information such as batch ID, QC values
  if(transcripts2genes==TRUE){
    if(biomart_cache){
      ensembl = biomaRt$useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset=biomart_dataset, host = "jul2015.archive.ensembl.org")
    } else {
      ensembl<-biomaRt$useMart("ensembl",dataset=biomart_dataset)
    }
    id_map<-biomaRt$getBM(c("ensembl_transcript_id","ensembl_gene_id"),mart=ensembl)
  } else {
    id_map<-NULL
  }
  #computationally expensive step- use parallel to speed up:
  ncores<-min(length(sample_list),parallel$detectCores())
  #d<-do.call("rbind",parallel$mclapply(sample_list,function(x){load_sample(x,id_map,new_cell_id=new_cell_id)},mc.cores=ncores))
  d<-do.call("rbind",parallel$mclapply(sample_list,load_sample,id_map,new_cell_id,tpm,mc.cores=ncores))
  vals<-if(tpm){ d$tpm }else{ d$est_counts } 
  m<-coo2matrix(d$gene_id,d$cell_id,vals)
  meta<-do.call("rbind",lapply(sample_list,qc2meta,new_cell_id))
  meta$pct_aligned<-with(meta,aligned_reads/total_reads)
  det<-data.frame(cell_id=colnames(m),detection=colSums(m>0))
  meta<-plyr$join(meta,det,by="cell_id")
  #sample_id_map<-sqldf("select distinct cell,cell_id,sample from d")
  mget(c("m","meta"))
}

# sample usage:
# tumors<-c("bc1","bc2","BC81-P6","BC-84-P12","BC-P2-0903")
# system.time(dat<-load_all_samples(tumors))
# m<-dat[["m"]]
# meta<-dat[["meta"]]
# all(colnames(m)%in%meta$cell_id)
# all(meta$cell_id%in%colnames(m))

#time: regular lapply: 118 seconds
#time: parallel mclapply with 5 cores: 55 seconds

