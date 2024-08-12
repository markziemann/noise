library("parallel")
library("edgeR")
library("DESeq2")
library("limma")
library("stringi")

########################################
# get some counts
########################################
countData<-function() {
# a is orig expression data
a<-read.table("https://raw.githubusercontent.com/markziemann/simde/master/start_data/ERR2539161/ERR2539161.se.tsv")

# GI is the gene information. It has gene name information for mapping accessions to gene symbols in GMT files
gi<-read.table("https://raw.githubusercontent.com/markziemann/simde/master/start_data/ERR2539161/GeneInfo.tsv",header=T,row.names=1)

# merge gene names
aa<-merge(a,gi,by=0)
aa<-aa[,-c(4:7)]
aa<-aggregate(. ~ GeneSymbol,aa,function(x) sum(as.numeric(as.character(x))))
aa$Row.names=NULL
rownames(aa)<-aa$GeneSymbol
aa$GeneSymbol=NULL
a<-aa[which(aa$ERR2539161>=10),,drop=F]
a
}

########################################
# generate some gene sets
########################################
randomGeneSets<-function(a,setsize,nsets){
gsets<-sapply( rep(setsize,nsets) , function(x) {list(as.character(sample(rownames(a),x))) } )
names(gsets)<-stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")
gsets
}

########################################
# simulate some gene expression data
########################################
simrna<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets) {

# N_REPS=5 ; SUM_COUNT=10000000 ; VARIANCE=0.2 ; FRAC_DE=0.05 ; FC=1 ; GMT="ReactomePathways.gmt"

library("edgeR")

df = NULL
for (k in paste0("data",1:(N_REPS*2)))  {
        b<-thinCounts(a,target.size=SUM_COUNT)
        colnames(b)=k
        df = cbind(df,b)
     }

# now need to only include gsets with 10 members in the 
gsets_sub<-which(unlist( lapply(gsets,function(x) { length(which(rownames(a) %in% as.character(unlist(x)))) >10 }  ) ) )
gsets<-gsets[which(names(gsets) %in% names(gsets_sub))]

#Number of differential gene sets
NDIF=round(length(gsets)*FRAC_DE)

if (VARIANCE>0) {
  #create some random values centred around 1 with some% error
  rand<-matrix(log2(rnorm(nrow(a)*N_REPS*2 , 2, VARIANCE)),ncol=N_REPS*2)
  #incorporate the noise
  df<-round(df*rand)
  #set any negative counts to zero
  df<-apply(df, 2, function(x) {ifelse(x < 0, 0, x)})
} 

if (NDIF>0) {
  message("prep fold changes")
  #Make even
  if ( NDIF%%2==1 ) { print("odd") ; NDIF=NDIF-1 }

  # sample some pathways to fiddle with
  DE_LIST<-sample(gsets , NDIF)

  # divide the list in 2 with half up and half down
  UP_LIST=sample(DE_LIST , NDIF/2)
  DN_LIST<-DE_LIST[!(DE_LIST %in% UP_LIST)]

  # now find a list of genes inside the pathways
  UP_DE<-unique(unlist(unname(UP_LIST)))

  # select the ones that are also in the profile
  UP_DE<-UP_DE[which(UP_DE %in% row.names(df))]
  #UP_DE<-sample(UP_DE,length(UP_DE)/2)

  # same for down genes
  DN_DE<-unique(unlist(unname(DN_LIST)))
  DN_DE<-DN_DE[which(DN_DE %in% row.names(df))]
  #DN_DE<-sample(DN_DE,length(DN_DE)/2)

  ITX<-intersect(UP_DE,DN_DE)
  # need to eliminate the overlapping ones for simplicity
  UP_DE<-setdiff(UP_DE,ITX)
  DN_DE<-setdiff(DN_DE,ITX)

  #reformat as df and add fold change
  UP_DE<-as.data.frame(UP_DE)
  UP_DE$V1<-2^FC
  colnames(UP_DE)=c("Gene","FC")
  DN_DE<-as.data.frame(DN_DE)
  DN_DE$V1<-2^-FC
  colnames(DN_DE)=c("Gene","FC")
  ALL_DE<-rbind(DN_DE,UP_DE)
  #Go back to list for downstream work
  UP_DE<-UP_DE$Gene
  DN_DE<-DN_DE$Gene
  NON_DE<-as.data.frame(setdiff(rownames(df),ALL_DE$Gene))
  colnames(NON_DE)="Gene"
  NON_DE$FC=1
  ALL_DE<-rbind(ALL_DE,NON_DE)
  ALL_DE<-ALL_DE[ order(as.vector(ALL_DE$Gene)) , ]
  message("incorporate changes")

  df <- df[ order(row.names(df)), ]
  df2<-cbind(df,ALL_DE)
  df2$Gene=NULL
} else {
  df2<-as.data.frame( df )
  df2$FC<- 1
  UP_DE=NULL
  DN_DE=NULL
  UP_LIST=NULL
  DN_LIST=NULL
}
ODD_COLS=(1:(ncol(df2)-1))[c(TRUE,FALSE)]
EVEN_COLS=(1:(ncol(df2)-1))[c(FALSE,TRUE)]
controls<-df2[,ODD_COLS]
colnames(controls)=paste0( "ctrl_" ,1:ncol(controls) )
treatments<-round(df2[,EVEN_COLS]*df2$FC)
colnames(treatments)=paste0( "trt_" ,1:ncol(treatments) )
x<-cbind(controls,treatments)
rownames(x)=rownames(df2)
#filter out genes that are not expressed
x<- x[which(rowSums(x)/ncol(x)>10),]
UP_DE<-intersect(UP_DE,rownames(x))
DN_DE<-intersect(DN_DE,rownames(x))
xx <- list("x" = x,"UP_DE"=UP_DE,"DN_DE"=DN_DE,"UP_LIST"=UP_LIST,"DN_LIST"=DN_LIST)
xx
}

#################################################
# a parallel repeat function
##################################################
#Thanks Gray Calhoun gcalhoun@iastate.edu for the following function
RepParallel <- function(n, expr, simplify = "array",...) {
      answer <-
        mclapply(integer(n), eval.parent(substitute(function(...) expr)),...)
      if (!identical(simplify, FALSE) && length(answer)) 
        return(simplify2array(answer, higher = (simplify == "array")))
      else return(answer)
    }
# RepParallel usage
#xxx<-RepParallel(10,simrna(a,5,10000000,0.2,20), simplify=F, mc.cores = detectCores() )


#################################################
# define DGE function
##################################################
dge_func <- function(x) {

  label = "simulate"
  y <- x[[1]]
  samplesheet <- as.data.frame(colnames(y))
  colnames(samplesheet) = "sample"
  samplesheet$trt <- factor(as.numeric(grepl("trt",colnames(y))))

  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData = y, colData = samplesheet, design = ~ trt )
  res <- DESeq(dds)
  z <- DESeq2::results(res)
  z <- as.data.frame(z[order(z$pvalue),])
  up <- rownames(subset(z,log2FoldChange>0 & padj<0.05))
  if (length(up)<200) { up <- rownames(head(subset(z,log2FoldChange>0 ),200)) }
  dn <- rownames(subset(z,log2FoldChange<0 & padj<0.05))
  if (length(dn)<200) { dn <- rownames(head(subset(z,log2FoldChange<0 ),200)) }
  x[[6]] <- z

  # edgeR
  group <- samplesheet$trt
  yy <- DGEList(counts=y,group=group)
  yyy <- normLibSizes(yy)
  design <- model.matrix(~group)
  fit <- glmQLFit(yyy,design)
  qlf <- glmQLFTest(fit,coef=2)
  x[[7]] <- topTags(qlf,number=Inf)

  # voom
  x
}

##################################
# aggregate function
##################################
agg_dge <- function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets, cores = 16) {

#TEST# N_REPS=5 ; SUM_COUNT=30000000 ; VARIANCE=0.45 ; FRAC_DE=0.05 ; FC=1 ; SIMS=8 ; DGE_FUNC="deseq2" ; gsets=gsets

xxx <- RepParallel(SIMS,simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets), simplify=F, mc.cores = cores )

# run deseq2
xxx <- mclapply(xxx , dge_func , mc.cores = cores )



# return the result
g=list()
for (f in 1:length(attributes(xxx))) {
  LABEL <- names(attributes(xxx)[f])
  LABEL <- as.data.frame(LABEL)
  g[[f]]<-cbind(unname(attributes(xxx)[f]),LABEL)
}
g
}
# x<-agg_dge(a,10,40000000,0.4,0.2,1,10,"deseq2",gsets) 

