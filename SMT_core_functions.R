library("gmp")
library("VennDiagram")
library("utils")
library('ape')
library('Rsamtools')
library('agrmt')

source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')

#filter microbial containing spots by UMI count
filterMicrobeSpot<-function(msrt,feature,min.count){
#  message(paste0("Searching ",feature))
  allfeaturedata=FetchData(msrt,vars=feature) %>% rowSums()
  targetspots=allfeaturedata[allfeaturedata>=min.count] %>% names()
  return(targetspots)
}

filterFeatures<-function(srt,min.occur=2,min.count=1,assay="Phyloseq_level"){
  DefaultAssay(srt)=assay
  data=FetchData(srt,cells=Cells(srt),vars=Features(srt))
  min.occur.counts=(data>=min.count) %>% colSums()
  validFeatures=names(min.occur.counts)[min.occur.counts>=min.occur]
  return(validFeatures)
}

CreateSeuratFromSMT.kraken<-function(SMT.out.path,taxa.level="genus",mode="combined"){
  #create seurat object
  raw=read.delim(paste0(SMT.out.path,'/matrix.mdx'),row.names = 1)
  colnames(raw)=colnames(raw) %>% gsub(".","-",fixed = T,.)
  tSeurat=Load10X_Spatial(SMT.out.path,filter.matrix=F,filename = "raw_feature_bc_matrix.h5")
  allbarcodes=Cells(tSeurat)
  invalidbarcodes=allbarcodes[!(allbarcodes %in% colnames(raw))]
  zeros=data.frame(row.names = row.names(raw),matrix(0,nrow = nrow(raw),ncol = length(invalidbarcodes)))
  colnames(zeros)=invalidbarcodes
  raw=cbind(raw,zeros)
  mSeurat=CreateSeuratObject(raw,assay="Spatial")
  mSeurat@images=tSeurat@images
  mSeurat@images$slice1@coordinates=mSeurat@images$slice1@coordinates[allbarcodes,]
  #create phyloseq object
  allfeatures=raw
  taxadb=read.delim(paste0(SMT.out.path,'/possorted_genome_bam.merged.microbe.tsv'),header=F)
  taxadb=taxadb[,c(8,9:16)]
  colnames(taxadb)<-c("taxid","superkingdom","kingdom","phylum","class","order","family","genus","species")
  taxadb[taxadb$kingdom=="kingdom_","kingdom"]<-taxadb[taxadb$kingdom=="kingdom_","superkingdom"] %>% gsub("super","",.)
  taxadb[taxadb$kingdom=="kingdom_","kingdom"]<-"kingdom_Unknown"
  taxadb[taxadb$phylum=="phylum_","phylum"]<-"phylum_Unknown"
  taxadb[taxadb$class=="class_","class"]<-"class_Unknown"
  taxadb[taxadb$order=="order_","order"]<-"order_Unknown"
  taxadb[taxadb$family=="family_","family"]<-"family_Unknown"
  taxadb[taxadb$genus=="genus_","genus"]<-"genus_Unknown"
  taxadb=taxadb %>% unique()
  rownames(taxadb)=taxadb$taxid
  taxadb=taxadb[,colnames(taxadb)!="taxid" & colnames(taxadb)!="superkingdom"]
  taxahier=taxadb[rownames(allfeatures),]
  phylo=phyloseq(otu_table(allfeatures,taxa_are_rows = T),tax_table(taxahier %>% as.matrix()))
  meta=data.frame(row.names = phylo@otu_table@.Data %>% colnames(),spot=phylo@otu_table@.Data %>% colnames())
  sample_data(phylo)<-meta
  phylo_level=taxa_level(phylo,taxa.level)
  phylo_count<-phylo_level@otu_table@.Data %>% t
  assay_phyloseq<-CreateAssayObject(phylo_count)
  if(mode=="microbe"){
    mSeurat@assays$Phyloseq_level=assay_phyloseq
    mSeurat@assays$Phyloseq_level@key="spatial_"
    mSeurat@assays$Phyloseq_raw=phylo
    DefaultAssay(mSeurat)="Phyloseq_level"
    return(mSeurat)
  }else if(mode=="combined"){
    cSeurat=Load10X_Spatial(SMT.out.path,filter.matrix=F,filename="raw_feature_bc_matrix.h5")
    cSeurat@assays$Phyloseq_level=assay_phyloseq
    cSeurat@assays$Phyloseq_level@key="spatial_"
    cSeurat@assays$Phyloseq_raw=phylo
    DefaultAssay(cSeurat)="Phyloseq_level"
    cSeurat=AddMetaData(cSeurat,metadata=colSums(cSeurat),col.name="nCount_Microbe")
    nfeature=colSums(cSeurat@assays$Phyloseq_level@counts>0)
    cSeurat=AddMetaData(cSeurat,metadata=nfeature,col.name="nFeature_Microbe")
    DefaultAssay(cSeurat)="Spatial"
    return(cSeurat)
  }else{
    stop("Invalid mode, microbe or combined")
  }
}

CreateSeuratFromSMT<-function(SMT.out.path,taxa.level="genus",mode="combined"){
  #create seurat object
  raw=read.delim(paste0(SMT.out.path,'/matrix.mdx'),row.names = 1)
  raw$barcode=NULL
  colnames(raw)=colnames(raw) %>% gsub(".","-",fixed = T,.)
  tSeurat=Load10X_Spatial(SMT.out.path,filter.matrix = F,filename = "raw_feature_bc_matrix.h5")
  allbarcodes=Cells(tSeurat)
  invalidbarcodes=allbarcodes[!(allbarcodes %in% colnames(raw))]
  zeros=data.frame(row.names = row.names(raw),matrix(0,nrow = nrow(raw),ncol = length(invalidbarcodes)))
  colnames(zeros)=invalidbarcodes
  raw<-cbind(raw,zeros)
  mSeurat=CreateSeuratObject(raw,assay="Spatial")
  mSeurat@images=tSeurat@images
  mSeurat@images$slice1@coordinates=mSeurat@images$slice1@coordinates[allbarcodes,]
  #create phyloseq object
  allfeatures=raw
  taxadb=read.delim(paste0(SMT.out.path,'/possorted_genome_bam.merged.microbe.tsv'),header=F)
  taxadb=taxadb[,c(4,(ncol(taxadb)-7):ncol(taxadb))]
  colnames(taxadb)<-c("taxid","superkingdom","kingdom","phylum","class","order","family","genus","species")
  taxadb[taxadb$kingdom=="kingdom_","kingdom"]<-taxadb[taxadb$kingdom=="kingdom_","superkingdom"] %>% gsub("super","",.)
  taxadb[taxadb$kingdom=="kingdom_","kingdom"]<-"kingdom_Unknown"
  taxadb[taxadb$phylum=="phylum_","phylum"]<-"phylum_Unknown"
  taxadb[taxadb$class=="class_","class"]<-"class_Unknown"
  taxadb[taxadb$order=="order_","order"]<-"order_Unknown"
  taxadb[taxadb$family=="family_","family"]<-"family_Unknown"
  taxadb[taxadb$genus=="genus_","genus"]<-"genus_Unknown"
  taxadb=taxadb %>% unique()
  rownames(taxadb)=taxadb$taxid
  taxadb=taxadb[,colnames(taxadb)!="taxid" & colnames(taxadb)!="superkingdom"]
  taxahier=taxadb[rownames(allfeatures),]
  phylo=phyloseq(otu_table(allfeatures,taxa_are_rows = T),tax_table(taxahier %>% as.matrix()))
  meta=data.frame(row.names = phylo@otu_table@.Data %>% colnames(),spot=phylo@otu_table@.Data %>% colnames())
  sample_data(phylo)<-meta
  phylo_level=taxa_level(phylo,taxa.level)
  phylo_count<-phylo_level@otu_table@.Data %>% t
  assay_phyloseq<-CreateAssayObject(phylo_count)
  if(mode=="microbe"){
    mSeurat@assays$Phyloseq_level=assay_phyloseq
    mSeurat@assays$Phyloseq_level@key="spatial_"
    mSeurat@assays$Phyloseq_raw=phylo
    DefaultAssay(mSeurat)="Phyloseq_level"
    return(mSeurat)
  }else if(mode=="combined"){
    cSeurat=Load10X_Spatial(SMT.out.path,filter.matrix=F,filename="raw_feature_bc_matrix.h5")
    cSeurat@assays$Phyloseq_level=assay_phyloseq
    cSeurat@assays$Phyloseq_level@key="spatial_"
    cSeurat@assays$Phyloseq_raw=phylo
    DefaultAssay(cSeurat)="Phyloseq_level"
    cSeurat=AddMetaData(cSeurat,metadata=colSums(cSeurat),col.name="nCount_Microbe")
    nfeature=colSums(cSeurat@assays$Phyloseq_level@counts>0)
    cSeurat=AddMetaData(cSeurat,metadata=nfeature,col.name="nFeature_Microbe")
    DefaultAssay(cSeurat)="Phyloseq_level"
    return(cSeurat)
  }else{
    stop("Invalid mode, microbe or combined")
  }
}

CreateSeuratFromSMT.MRCA<-function(SMT.out.path,taxa.level="genus",mode="combined"){
  #create seurat object
  raw=read.delim(paste0(SMT.out.path,'/matrix.mdx'),row.names = 1)
  raw$barcode=NULL
  colnames(raw)=colnames(raw) %>% gsub(".","-",fixed = T,.)
  tSeurat=Load10X_Spatial(SMT.out.path,filter.matrix = F,filename = "raw_feature_bc_matrix.h5")
  allbarcodes=Cells(tSeurat)
  invalidbarcodes=allbarcodes[!(allbarcodes %in% colnames(raw))]
  zeros=data.frame(row.names = row.names(raw),matrix(0,nrow = nrow(raw),ncol = length(invalidbarcodes)))
  colnames(zeros)=invalidbarcodes
  raw<-cbind(raw,zeros)
  mSeurat=CreateSeuratObject(raw,assay="Spatial")
  mSeurat@images=tSeurat@images
  mSeurat@images$slice1@coordinates=mSeurat@images$slice1@coordinates[allbarcodes,]
  #create phyloseq object
  allfeatures=raw
  taxadb=read.delim(paste0(SMT.out.path,'/taxa.tsv'),header=F)
  colnames(taxadb)<-c("taxid","superkingdom","kingdom","phylum","class","order","family","genus","species")
  taxadb[taxadb$kingdom=="kingdom_","kingdom"]<-taxadb[taxadb$kingdom=="kingdom_","superkingdom"] %>% gsub("super","",.)
  taxadb[taxadb$kingdom=="kingdom_","kingdom"]<-"kingdom_Unknown"
  taxadb[taxadb$phylum=="phylum_","phylum"]<-"phylum_Unknown"
  taxadb[taxadb$class=="class_","class"]<-"class_Unknown"
  taxadb[taxadb$order=="order_","order"]<-"order_Unknown"
  taxadb[taxadb$family=="family_","family"]<-"family_Unknown"
  taxadb[taxadb$genus=="genus_","genus"]<-"genus_Unknown"
  taxadb=taxadb %>% unique()
  rownames(taxadb)=taxadb$taxid
  taxadb=taxadb[,colnames(taxadb)!="taxid" & colnames(taxadb)!="superkingdom"]
  taxahier=taxadb[rownames(allfeatures),]
  phylo=phyloseq(otu_table(allfeatures,taxa_are_rows = T),tax_table(taxahier %>% as.matrix()))
  meta=data.frame(row.names = phylo@otu_table@.Data %>% colnames(),spot=phylo@otu_table@.Data %>% colnames())
  sample_data(phylo)<-meta
  phylo_level=taxa_level(phylo,taxa.level)
  phylo_count<-phylo_level@otu_table@.Data %>% t
  assay_phyloseq<-CreateAssayObject(phylo_count)
  if(mode=="microbe"){
    mSeurat@assays$Phyloseq_level=assay_phyloseq
    mSeurat@assays$Phyloseq_level@key="spatial_"
    mSeurat@assays$Phyloseq_raw=phylo
    DefaultAssay(mSeurat)="Phyloseq_level"
    return(mSeurat)
  }else if(mode=="combined"){
    cSeurat=Load10X_Spatial(SMT.out.path,filter.matrix=F,filename="raw_feature_bc_matrix.h5")
    cSeurat@assays$Phyloseq_level=assay_phyloseq
    cSeurat@assays$Phyloseq_level@key="spatial_"
    cSeurat@assays$Phyloseq_raw=phylo
    DefaultAssay(cSeurat)="Phyloseq_level"
    cSeurat=AddMetaData(cSeurat,metadata=colSums(cSeurat),col.name="nCount_Microbe")
    nfeature=colSums(cSeurat@assays$Phyloseq_level@counts>0)
    cSeurat=AddMetaData(cSeurat,metadata=nfeature,col.name="nFeature_Microbe")
    DefaultAssay(cSeurat)="Phyloseq_level"
    return(cSeurat)
  }else{
    stop("Invalid mode, microbe or combined")
  }
}


SMTPlot<-function(object,features,images=NULL,crop=TRUE,slot="data",
      min.cutoff=NA,max.cutoff=NA,ncol=NULL,combine=TRUE,
      pt.size.factor=1.3,alpha=c(1, 1),stroke=0.25,interactive=FALSE,legend.title=NULL,
      information=NULL,colors=c("gray","red"),legend.size=4,legend.direction="horizontal",
      legend.position="top",legend.text.size=10,plot.margin=0.5,legend.margin=0.5){
  SpatialFeaturePlot(object=object,features=features,images=images,crop=crop,slot=slot,
      min.cutoff=min.cutoff,max.cutoff=max.cutoff,ncol=ncol,combine=combine,
      pt.size.factor=pt.size.factor,alpha=alpha,stroke=stroke,interactive=interactive,
      information=information)+scale_fill_gradientn(colours=colors,guide="legend")+
    theme(legend.direction=legend.direction,
          legend.position=legend.position,
          legend.margin=margin(legend.margin,legend.margin,legend.margin,legend.margin,"cm"),
          legend.text=element_text(size=legend.text.size),legend.title = element_text(size=legend.text.size),
          plot.margin=margin(plot.margin,plot.margin,plot.margin,plot.margin,"cm"))+
    guides(fill=guide_legend(title=legend.title,override.aes=list(size=legend.size)))
}

CallTaxaTableFromTSV<-function(mergedblastout.path,taxa.level,sample.id,scale.factor=100000){
  res=read.delim(mergedblastout.path,header=F)
  res=res %>% dplyr::filter(.,(V9 %in% c("superkingdom_Archaea","superkingdom_Bacteria","superkingdom_Viruses")) |
                              (V10 %in% c("kingdom_Fungi")))
  colnames(res)<-c("read","taxid1","target_seq","length","evalue","identity","coverage","taxid2","Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus",
  "Species")
  res=res[,taxa.level] %>% table
  res=data.frame(count=res %>% as.numeric(),sample=rep(sample.id,length(res)),taxa=names(res))
  if(scale.factor==0){
    res$count=res$count
  }else{
    res$count=res$count*scale.factor/sum(res$count)
  }
  return(res)
}

SMTlevel<-function(object,new.taxa.level){
  phylo=object@assays$Phyloseq_raw
  new_phylo_level=taxa_level(phylo,new.taxa.level)
  phylo_count=new_phylo_level@otu_table@.Data %>% t
  assay_phyloseq=CreateAssayObject(phylo_count)
  object@assays$Phyloseq_level=assay_phyloseq
  DefaultAssay(object)="Phyloseq_level"
  nfeature=colSums(object@assays$Phyloseq_level@counts>0)
  object=AddMetaData(object,metadata=nfeature,col.name="nFeature_Microbe")
  return(object)
}

SMTStackPlot<-function(object,taxa.threshold=0.1,min.nCount=10){
  DefaultAssay(object)="Phyloseq_level"
  df=FetchData(object,vars = Features(object)) %>% t %>% as.data.frame()
  validspot=colnames(df)[colSums(df)>=min.nCount]
  if(length(validspot)==0){
    stop("No spot meet the min nCount threshold.")
  }
  df=df[,validspot]
  mostEnrichedTaxa=rowSums(df) %>% sort() %>% tail(.,1) %>% names()
  order=(df[mostEnrichedTaxa,validspot]/colSums(df[,validspot])) %>% sort(.,decreasing=T) %>% names()
  df=TransformDataframe(df,taxa.threshold)
#  return(df)
  df=dplyr::filter(df,spot %in% validspot)
  df$spot=factor(df$spot,levels=order)
  taxaorder=(arrange(df[df$features!="others",],by=desc(values)))$features %>% append(mostEnrichedTaxa,.) %>% unique()
  taxaorder=append(taxaorder,"others")
  df$features=factor(df$features,levels=taxaorder)
#  return(taxaorder)
  message("plotting")
  ggplot(df)+geom_bar(aes(x=spot,y=values,fill=features),stat="identity",position="fill")+
    scale_fill_d3("category20")+theme_minimal()+theme(axis.text.x=element_blank())+
    xlab("")
}

TransformRawDf<-function(df,cutoff=0){
  data=df
  cns=colnames(data)
  rns=rownames(data)
  ret=data.frame("feat1"=character(0),"feat2"=character(0),"values"=numeric(0))
  for (i in cns){
    print(i)
    if (rns[data[,i] %>% abs() >= cutoff] %>% length() == 0){
      next
    }
    feat2=rns[data[,i] %>% abs() >= cutoff]
    values=data[feat2,i] %>% as.numeric()
    feat1=rep(i,length(feat2))
    subdf=data.frame(feat1=feat1,feat2=feat2,values=values)
    ret=rbind(ret,subdf)
  }
  return(ret)
}

#randomly putting a number of UMI into spots for permutation, calculate possibility of each permutation.
RandomUMIPlacing<-function(srt,nCount,n.permutation=1000){
  barcodes=Cells(srt)
#  out=NULL
  clnum=detectCores()/2
#  message(paste0(clnum," cores detected, will use ",clnum))
  cl=makeCluster(getOption("cl.cores", clnum))
  res=parLapply(cl,1:n.permutation,SingleRandomUMIPlacing,barcodes=barcodes,nCount=nCount)
  stopCluster(cl)
#  for(i in 1:n.permutation){
#    sample.result=sample(barcodes,nCount,T)
#    unique.barcode.num=length(table(sample.result))
#    possibility=1/(choose(4992,unique.barcode.num)*choose(nCount-1,unique.barcode.num-1))
#    way.placing.part1=prodmod(4992)/(prodmod(unique.barcode.num)*prodmod(4992-unique.barcode.num))
#    way.placing.part2=prodmod(nCount-1)/(prodmod(unique.barcode.num-1)*prodmod(nCount-unique.barcode.num))
#    way.placing=way.placing.part1*way.placing.part2
#    out=append(out,way.placing)
#  }
  return(res)
}

RandomSpotPlacing<-function(srt,nSpot,n.permutation=1000){
  source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')
  barcodes=Cells(srt)
  coords=srt@images$slice1@coordinates
  clnum=detectCores()/2
#  message(paste0(clnum," cores detected, will use ",clnum))
  cl=makeCluster(getOption("cl.cores", clnum))
  res=parLapply(cl,1:n.permutation,SingleRandomSpotSelecting,barcodes=barcodes,nSpot=nSpot,coords=coords)
  stopCluster(cl)
  res=res %>% unlist()
  return(res)
}

SingleRandomUMIPlacing<-function(index,barcodes,nCount){
  source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')
  sample.result=sample(barcodes,nCount,T)
  unique.barcode.num=length(table(sample.result))
  way.placing.part1=prodmod(4992)/(prodmod(unique.barcode.num)*prodmod(4992-unique.barcode.num))
  way.placing.part2=prodmod(nCount-1)/(prodmod(unique.barcode.num-1)*prodmod(nCount-unique.barcode.num))
  way.placing=way.placing.part1*way.placing.part2
  return(way.placing)
}

SingleRandomSpotSelecting<-function(index,barcodes,nSpot,coords){
#  source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')
  sample.result=sample(barcodes,nSpot,F)
  sample.coords=coords[sample.result,c("imagerow","imagecol")]
  res=median(sort(as.vector(dist(sample.coords))))
#  debug=list()
#  debug[[1]]=res
#  debug[[2]]=sample.result
#  debug[[3]]=sample.coords
  return(res)
#  return(debug)
}

CalSpotSignificance<-function(srt,n.permutation){
  pb=txtProgressBar(style=3)
  microbes=rownames(srt@assays$Phyloseq_level@counts)
  coords=srt@images$slice1@coordinates
  out=data.frame(taxa=NULL,nSpot=NULL,p.value=NULL)
  samples=list()
  i=1
  for(feature in microbes){
    barcodes=filterMicrobeSpot(srt,feature,1)
    nSpot=length(barcodes)
    if(nSpot==1){
      row=data.frame(taxa=feature,nSpot=1,p.value=NA)
    }else if(nSpot %in% out$nSpot){
      feature.median.dist=median(dist(coords[barcodes,c("imagerow","imagecol")]))
      sampled.median.dist=samples[[nSpot]]
      p.value=sum(sampled.median.dist<feature.median.dist)/n.permutation
      row=data.frame(taxa=feature,nSpot=nSpot,p.value=p.value)
    }else{
      feature.median.dist=median(dist(coords[barcodes,c("imagerow","imagecol")]))
      sampled.median.dist=RandomSpotPlacing(srt=srt,nSpot=nSpot,n.permutation=n.permutation)
      p.value=sum(sampled.median.dist<feature.median.dist)/n.permutation
      row=data.frame(taxa=feature,nSpot=nSpot,p.value=p.value)
      samples[[nSpot]]=sampled.median.dist
    }
    out=rbind(out,row)
    setTxtProgressBar(pb,i/length(microbes))
    i=i+1
  }
  return(out)
}

#calculate significance of microbe signals to distinguish contamination from true signals.
CalUMISignificance<-function(srt,n.permutation){
  pb=txtProgressBar(style=3)
  counts=srt@assays$Phyloseq_level@counts
  features=rownames(counts)
  out=data.frame(taxa=NULL,nCount=NULL,p.value=NULL)
  i=1
  sampled=list()
  for(feature in features){
#    message(paste0("Progress: ",i,"/",length(features)))
    barcode.unique=sum(counts[feature,]>0)
    umi=sum(counts[feature,])
#    possibility=1/(choose(4992,barcode.unique)*choose(umi-1,barcode.unique-1))
    way.placing.part1=prodmod(4992)/(prodmod(barcode.unique)*prodmod(4992-barcode.unique))
    way.placing.part2=prodmod(umi-1)/(prodmod(barcode.unique-1)*prodmod(umi-barcode.unique))
    way.placing=way.placing.part1*way.placing.part2
    if(length(sampled)<umi||(length(sampled)>=umi&is.null(sampled[[umi]]))){
      permutated.ways=RandomUMIPlacing(srt,nCount=umi,n.permutation=n.permutation)
      sampled[[umi]]=permutated.ways
    }else{
      permutated.ways=sampled[[umi]]
    }
    count=0
    for(n in 1:length(permutated.ways)){
      if(way.placing>=permutated.ways[[n]]){
        count=count+1
      }
    }
    p.value=count/n.permutation
    res=data.frame(taxa=feature,nCount=umi,unique.barcode.num=barcode.unique,p.value=p.value)
    out=rbind(out,res)
    setTxtProgressBar(pb,i/length(features))
    i=i+1
  }
  return(out)
}

prodmod<-function(x){
  if(x==0|x==1){
    y=1
  }else{
    y=prod(1:x %>% as.bigz())
  }
  return(y)
}

SMTPlotMod<-function(srt,feature,pt.size=1.1,colors=c("#FFCC00","red"),alpha=c(0.1,6),legend.title.override=NULL,
                     min.count=NULL,max.count=NULL,legend.text.size=13,log.scale=T,legend.position="top",
                     legend.direction="horizontal",legend.pt.size=3,on.tissue.only=F){
  coord=srt@images$slice1@coordinates
  if(on.tissue.only){
    count=FetchData(srt,vars=feature,cells=filterTissueSpots(srt)) %>% as.data.frame()
  }else{
    count=FetchData(srt,vars=feature) %>% as.data.frame()
  }
  colnames(count)="count"
  if(is.null(max.count)){
    max.count=max(count$count)
  }
  if(is.null(min.count)){
    min.count=min(count$count)
  }
  count[count$count<min.count,"count"]=min.count
  count[count$count>max.count,"count"]=max.count
  legend.title=feature %>% gsub("^.*-","",.)
  if(log.scale){
    count=log2(count+1)
    legend.title=paste0("log ",legend.title)
  }
  if(!is.null(legend.title.override)){
    legend.title=legend.title.override
  }
#  return(count)
  data=mergeByRownames(coord,count)
  ggplot()+geom_point(data=data %>% dplyr::filter(.,tissue==0),aes(x=imagecol,y=imagerow),alpha=0.1,size=pt.size)+
    geom_point(data=data %>% dplyr::filter(.,tissue==1),aes(x=imagecol,y=imagerow),alpha=0.3,size=pt.size)+
    geom_point(data=data %>% dplyr::filter(.,count>min.count),aes(x=imagecol,y=imagerow,alpha=count,color=count),size=pt.size)+scale_alpha_continuous(range=alpha)+
    scale_y_reverse()+coord_fixed()+scale_color_gradientn(colours = colors)+
    theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_blank(),
          plot.margin=margin(0.5,0.5,0.5,0.5,unit="cm"),
          legend.key=element_rect(fill = "white"),
          legend.text=element_text(size=legend.text.size),legend.title=element_text(size=legend.text.size),
          legend.margin=margin(0,0,0,0,unit="cm"),
          legend.position = legend.position,
          legend.direction = legend.direction)+
    guides(color=guide_legend(title=legend.title,
           override.aes=list(size=legend.pt.size)),alpha=guide_legend(title=legend.title))
}

denoise<-function(srt,n.permutation=1000,strict.mode=T,pvalue=0.01,return.raw=F){
  message("denoise by spot..")
#  valid1=CalSpotSignificance(srt,n.permutation)
  valid0=CalSpotSignificance(srt,n.permutation)
#  message("denoise by UMI..")
#  valid2=CalUMISignificance(srt,1000)
  valid1=valid0 %>% dplyr::filter(.,p.value<pvalue)
#  valid2=valid2 %>% dplyr::filter(.,p.value<pvalue)
  if(strict.mode){
#    valid.taxa=intersect(valid1$taxa,valid2$taxa)
    valid.taxa=valid1$taxa
  }else{
    valid.taxa=valid1$taxa
#    valid.taxa=unique(valid1$taxa,valid2$taxa)
  }
  if("Spatial" %in% names(srt@assays)){
    host.assay=GetAssay(srt,assay = "Spatial")
    srt=subset(srt,features=valid.taxa)
    srt@assays$Spatial=host.assay
  }else{
    srt=subset(srt,features=valid.taxa)
  }
  if(return.raw){
    return(valid0)
  }else{
    return(srt)
  }
}

Features<-function(srt){
  features=rownames(srt)
  return(features)
}

#smooth count of spots by averaging it with its perispots expression
#srt: Seurat Object as input
#weight.cent: how much proportion of the spot's expression should be add to the sum
#weight.peri: how much proportion of perispots' expression should be add to the sum
#mode: "default", "smt"
#tech: only support "visium"
SmoothSingleSpot<-function(srt,spot,weight.cent=1,weight.peri=0.5,mode="default",tech="visium"){
  if(tech=="visium"){
    if(mode=="default"){
      expression.old=srt@assays$Spatial@counts
    }else if(mode=="smt"){
      expression.old=srt@assays$Phyloseq_level@counts
    }else{
      stop("Mode you specified is not supported, it should be 'smt' or 'default'.")
    }
    peri.spot=SelectNeighbors(srt=srt,spot=spot,layout="hex")
    expression.new=(expression.old[,spot]*weight.cent+rowSums(expression.old[,peri.spot]*weight.peri))/(1*weight.cent+6*weight.peri)
    expression.new=as.data.frame(expression.new)
    colnames(expression.new)=spot
    return(expression.new)
  }else{
    stop("Other tech not supported")
  }
}

SmoothSpots<-function(srt,weight.cent=1,weight.peri=0.5,mode="default",tech="visium"){
  message("Smoothing with lapply ...")
  res=lapply(FUN=SmoothSingleSpot,X=Cells(srt),srt=srt,weight.cent=weight.cent,weight.peri=weight.peri,mode=mode,tech=tech)
  message("Merging results into single dataframe ...")
  res=res %>% as.data.frame()
  colnames(res)=colnames(res) %>% gsub(".","-",.,fixed=T)
  message("Creating new assay ...")
  srt@assays$Smooth=CreateAssayObject(counts=res)
  srt=AddMetaData(srt,metadata=srt@assays$Smooth@counts %>% colSums(),col.name="nCount_Smooth")
  srt=AddMetaData(srt,metadata=(srt@assays$Smooth@counts>0) %>% colSums(),col.name="nFeature_Smooth")
  return(srt)
}

denoiseByMoran.I.test<-function(srt,feature){
  count=FetchData(srt,vars=feature)
  coords=srt@images$slice1@coordinates[,c("imagerow","imagecol")]
  data=mergeByRownames(coords,count)
  return(data)
  dists=as.matrix(dist(data[,c("imagerow","imagecol")]))
  dists.inv=1/dists
  diag(dists.inv)=0
  dists.inv[is.infinite(dists.inv)]=0
  re=Moran.I(data[,feature],dists.inv,alternative='g')
  return(re)
}

findSingleSMTcor<-function(msrt,microbe,max.layer=2,min.gene.occurence=5,min.gene.umi.count=50){
  DefaultAssay(msrt)="Phyloseq_level"
  spots=filterMicrobeSpot(msrt,microbe,min.count=1)
  if(length(spots)<5){
    stop("There's not enough spots for correlation for this taxa")
  }
  expressions=list()
  validmicrobespots=NULL
  for(spot in spots){
    perispots=SelectNeighborLayers(srt=msrt,spot=spot,layer_num=max.layer,layout="hex")
    tissuespots=filterTissueSpots(msrt,spots=c(perispots,spot))
    if(length(tissuespots)==0){
      next
    }
    DefaultAssay(msrt)="Spatial"
#    geneexpressions=colSums(FetchData(msrt,vars=Features(msrt),cells=tissuespots))/length(tissuespots)
    geneexpressions=FetchData(msrt,vars=Features(msrt),cells=tissuespots) %>% map(.,median) %>% unlist()
    expressions[[spot]]=geneexpressions
    validmicrobespots=append(validmicrobespots,spot)
    DefaultAssay(msrt)="Phyloseq_level"
  }
  re=base::Reduce(f=rbind,x=expressions)
  rownames(re)=validmicrobespots
  re=re[,colSums(re>0)>=min.gene.occurence]
  re=re[,colSums(re)>=min.gene.umi.count]
  if(ncol(re)==0){
    return(NULL)
  }
  richness=FetchData(msrt,vars=microbe,cells=validmicrobespots)
  re=mergeByRownames(richness,re)
  re=log2(re+1)
#  return(re)
  mat1=re[,1] %>% as.data.frame()
  mat2=re[,2:ncol(re)] %>% as.data.frame()
  cors=matcor(mat1,mat2)
  return(cors)
}

findSMTcor<-function(msrt,features,max.layer=2,min.gene.occurence=5,min.gene.umi.count=50){
  outs=list()
  for(feature in features){
    if(length(filterMicrobeSpot(msrt,feature,1))<5){
      message(paste0("Feature ",feature," did not meet limit to perform correlation analysis"))
      next
    }
    out=findSingleSMTcor(msrt,feature,max.layer=max.layer,min.gene.occurence=min.gene.occurence,min.gene.umi.count=min.gene.umi.count)
    if(is.null(out)){
      message(paste0("Feature ",feature," did not meet limit to perform correlation analysis"))
      next
    }
    outs[[feature]]=out
  }
  return(outs)
}

#divide the tissue into different part if multi-tissue were laid on this area, mode can be "SMT" to obtain peri-partition.
getPartition<-function(srt,partition.num="auto",eps=2,MinPts=3,labels="auto",mode="auto",SMT.layer=2){
  data=srt@images$slice1@coordinates
  data=dplyr::filter(data,tissue==1)
  out=dbscan(data[,c("row","col")],eps=eps,MinPts=MinPts)
  names(out$cluster)=rownames(data)
  ressrt=AddMetaData(srt,metadata=out$cluster,col.name="partition")
  if(is.numeric(partition.num)){
    target.cluster=table(out$cluster) %>% sort() %>% tail(.,partition.num) %>% names()
    ressrt@meta.data[!(ressrt@meta.data$partition %in% target.cluster),"partition"]="non-tissue spots"
  }
  if(is.vector(labels)&length(labels)>=2){
    if(is.numeric(partition.num)&length(labels)==partition.num){
      oldlabel=table(ressrt$partition) %>% names()
      i=1
      while(i<length(oldlabel)){
        ressrt@meta.data[ressrt@meta.data$partition==oldlabel[i],"partition"]=labels[i]
        i=i+1
      }
      newlevel=append(labels,"non-tissue spots")
      ressrt@meta.data$partition=factor(ressrt@meta.data$partition,levels=newlevel)
    }else{
      message("The number of labels you have provided is not the same as partition number you specified")
      stop()
    }
  }
  if(mode=="SMT"){
    #    return(ressrt)
    oldlabels=table(ressrt@meta.data$partition) %>% names()
    oldlabels=oldlabels[oldlabels!="non-tissue spots"]
    peritargets=list()
    i=1
    for (label in oldlabels){
      target=Cells(subset(ressrt,partition==label))
      peritarget=SelectNeighborLayers(srt=ressrt,spots=target,layer_num=SMT.layer)
      peritargets[[i]]=peritarget
      i=i+1
      ressrt@meta.data[unique(append(target,peritarget)),"SMTpartition"]=label
    }
    sharedspots=unique(unlist(peritargets)[duplicated(unlist(peritargets))])
    ressrt@meta.data[!(ressrt@meta.data$SMTpartition %in% oldlabels),"SMTpartition"]="non-SMT spots"
    ressrt@meta.data[sharedspots,"SMTpartition"]="ambiguous region"
  }
  return(ressrt)
}

#for multi-partiton sample, subset by partition, each subset can be defined as rectangular or partition shaped region.
subsetPartition<-function(srt,part,margin=2,mode="rec"){
  if(!("partition" %in% colnames(srt@meta.data))){
    stop("partition info not found in metadata, run getPartition first")
  }
  partitionspots=Cells(subset(srt,partition==part))
  coord=srt@images$slice1@coordinates[partitionspots,]
  if(mode=="rec"){
    max.row=max(coord$row)+1*margin
    min.row=min(coord$row)-1*margin
    max.col=max(coord$col)+2*margin
    min.col=min(coord$col)-2*margin
    subsetspots=dplyr::filter(srt@images$slice1@coordinates,row<=max.row,row>=min.row,col<=max.col,col>=min.col) %>% rownames()
  }else if(mode=="shape"){
    subsetspots=SelectNeighborLayers(srt,partitionspots,layer_num=margin)
    subsetspots=append(subsetspots,partitionspots) %>% unique()
  }
  ressrt=subset(srt,cells=subsetspots)
  #SpatialDimPlot(ressrt,"partition")
  return(ressrt)
}

filterFeatures<-function(srt,min.count=0,min.spot.occurrence=3){
  spot.occurence=FetchData(srt,vars=Features(srt))
  spot.occurence=colSums(spot.occurence>=min.count)
  spot.occurence=spot.occurence[spot.occurence>=min.spot.occurrence]
  out=names(spot.occurence)
  return(out)
}

findPairedFeatureCor<-function(srt,feature1,feature2,min.coocurrence=4,use.Smooth=F){
  if(!use.Smooth){
    f1spots=filterMicrobeSpot(srt,feature=feature1,min.count=1)
    f2spots=filterMicrobeSpot(srt,feature=feature2,min.count=1)
    if(length(intersect(f1spots,f2spots))<min.coocurrence){
      message(paste0(feature1," and ",feature2," do not have enough shared spots for correlation analysis, return 0"))
      out=data.frame(feat1=feature1,feat2=feature2,values=0)
      return(out)
    }else{
      data=FetchData(srt,vars=c(feature1,feature2),cells=intersect(f1spots,f2spots))
      out=as.data.frame(cor(data))
      out=TransformRawDf(out[1],cutoff=0)
      out=out[2,]
      return(out)
    }
  }
}

findOnevsNFeatureCor<-function(srt,feature,features,min.coocurrence=4,use.Smooth=F){
  outs=data.frame(feat1=character(0),feat2=character(0),values=numeric(0))
  for (feat in features) {
    out=findPairedFeatureCor(srt,feature,feat,min.coocurrence=min.coocurrence,use.Smooth=use.Smooth)
    if(!is.na(out$feat2)){
      outs=rbind(outs,out)
    }
  }
  rownames(outs)=1:nrow(outs)
  return(outs)
}

findNvsNFeatureCor<-function(srt,features,min.occur=2,min.count=1,font.size=1,
                             assay="Phyloseq_level",mask="genus-",return.matrix=F){
  DefaultAssay(srt)=assay
  validfeatures=filterFeatures(srt,min.spot.occurrence=min.occur,min.count=min.count)
  wantedvalidfeatures=intersect(features,validfeatures)
  wantedinvalidfeatures=setdiff(features,wantedvalidfeatures)
  message(paste0("These features did not meet the cut-off set or doesn't exist: ",paste(wantedinvalidfeatures,collapse=",")))
  data=FetchData(srt,vars=wantedvalidfeatures,cells=filterMicrobeSpot(srt,feature=wantedvalidfeatures,min.count=1))
  cors=cor(data)
  rownames(cors)<-rownames(cors) %>% gsub(mask,"",.)
  colnames(cors)<-colnames(cors) %>% gsub(mask,"",.)
  if(return.matrix){
    return(cors)
  }
  corrplot(cors,method = "circle",order="hclust",type="lower",tl.col = "black",
      cl.cex = font.size,tl.cex = font.size,tl.srt=45,col=colorRampPalette(c("blue", "white", "red"))(100))
}

igvplot<-function(sam.path,ref.path,thickness=0.15,ymin=0,space=0.05,read.color="darkgrey",
                  ref.color="black",view.part=NULL,show.entire.ref=F,return.snp=F,return.reads.info=F,
                  return.track.info=F){
  ref=read_lines(ref.path)
  ref=ref[2:length(ref)] %>% paste(.,collapse = "")
  ref.length=nchar(ref)
  self.sam=read.delim(sam.path,header = F)
  if(is.null(view.part)){
    self.sam=self.sam[,c(1,3,4,6,10)]
  }else{
    self.sam=self.sam[view.part,c(1,3,4,6,10)]
  }
  colnames(self.sam)=c("read","ref","start","CIGAR","seq")
#  self.sam$extended=self.sam$CIGAR %>% gsub("[[:digit:]]*[SH]","",.) %>%
#  strsplit(.,"[MDI]") %>% lapply(.,as.numeric) %>% lapply(.,sum) %>% unlist()
  self.sam$extended=self.sam$CIGAR %>% gsub("[[:digit:]]*[SHI]","",.) %>%
    strsplit(.,"[MD]") %>% lapply(.,as.numeric) %>% lapply(.,sum) %>% unlist()
  self.sam$end=self.sam$start+self.sam$extended-1
  self.sam=self.sam[,c("read","ref","start","end","extended","CIGAR","seq")]
  self.sam=self.sam %>% arrange(.,by=start)
  if(return.reads.info){
    return(self.sam)
  }
  self.sam[,c("ymin","ymax","track")]=NA
  self.sam[1,"ymin"]=ymin
  self.sam[1,"ymax"]=ymin+thickness
  self.sam[1,"track"]=1
  message("generating track data...")
  for(i in 2:nrow(self.sam)){
    suitable.track=1
#    message(i)
#    if(i==2547){
#      return(self.sam)
#    }
#    message(max(self.sam[self.sam$track==suitable.track,"end"]
#                [!is.na(self.sam[self.sam$track==suitable.track,"end"])]))
#    message(self.sam$start[i])
    while(suitable.track<=max(self.sam$track[1:(i-1)])&
          max(self.sam[self.sam$track==suitable.track,"end"]
              [!is.na(self.sam[self.sam$track==suitable.track,"end"])])>=self.sam$start[i]){
      suitable.track=suitable.track+1
      if(suitable.track>max(self.sam$track[1:(i-1)])){
        break
      }
    }
    self.sam$ymin[i]=ymin+(suitable.track-1)*(thickness+space)
    self.sam$ymax[i]=self.sam$ymin[i]+thickness
    self.sam$track[i]=suitable.track
  }
  #form indel plot data
  message("forming indel plot data...")
  indel.data=self.sam[,c("CIGAR","start","end","ymin","ymax")]
  indel.data=indel.data[grep(pattern="[DI]",x=indel.data$CIGAR),]
  indel.plotdata=data.frame("CIGAR"=character(0),"start"=numeric(0),"ymin"=numeric(0),
                            "ymax"=numeric(0),"end"=numeric(0),"indel"=character(0))
  if(nrow(indel.data)!=0){
    for(i in 1:nrow(indel.data)){
      indel.num=length(unlist(strsplit(x=indel.data$CIGAR[i],"[DI]")))-1
      indel.plotdata.sub=data.frame(
        CIGAR=rep(indel.data$CIGAR[i],indel.num),
        start=rep(indel.data$start[i],indel.num),
        ymin=rep(indel.data$ymin[i],indel.num),
        ymax=rep(indel.data$ymax[i],indel.num))
      valid.extension=indel.data$CIGAR[i] %>% gsub("[[:digit:]]*[SH]","",.)
      valid.extension=valid.extension %>% gsub("[[:digit:]]*I","0I",.)
      units=unlist(str_extract_all(valid.extension,"[[:digit:]]*[MDI]"))
      indel.index=grep("[ID]",units)
      units=gsub("[MDI]","",units) %>% as.numeric()
      extends=cumsum(units)[indel.index-1]
      indel.plotdata.sub$start=indel.plotdata.sub$start+extends
      indel.plotdata.sub$end=indel.plotdata.sub$start
      indel.plotdata.sub$indel=unlist(str_extract_all(valid.extension,"[DI]"))
      indel.plotdata.sub[indel.plotdata.sub$indel=="D","indel"]="deletion"
      indel.plotdata.sub[indel.plotdata.sub$indel=="I","indel"]="insertion"
      indel.plotdata=rbind(indel.plotdata,indel.plotdata.sub)
    }
  }

  #form SNP plot data
  message("forming SNP data...")
  snp.plotdata=call.snp.from.sam(self.sam,ref)

  #form ref plot data
  ref.row=data.frame("read"="ref","ref"="ref","start"=min(self.sam$start),"end"=max(self.sam$end),"extended"=max(self.sam$end),
                      "CIGAR"="ref","seq"="ref","ymin"=ymin-0.2-1.5*thickness,"ymax"=ymin-0.2)
  if(show.entire.ref){
    ref.row=data.frame("read"="ref","ref"="ref","start"=1,"end"=ref.length,"extended"=ref.length,
                       "CIGAR"="ref","seq"="ref","ymin"=ymin-0.2-1.5*thickness,"ymax"=ymin-0.2)
  }

  #plot
  p=ggplot()+geom_rect(data=self.sam,aes(xmin=start-1,xmax=end,ymin=ymin,ymax=ymax),fill=read.color,alpha=0.5)+
    geom_rect(data=ref.row,aes(xmin=start-1,xmax=end,ymin=ymin,ymax=ymax),fill=ref.color)+
    geom_rect(data=indel.plotdata,aes(xmin=start-1,xmax=end,ymin=ymin,ymax=ymax,fill=indel))+
    geom_rect(data=snp.plotdata,aes(xmin=start-1,xmax=end,ymin=ymin,ymax=ymax,fill=snp))+
    ylim(c(ref.row$ymin-5,0.5+max(self.sam$ymax)))+
    scale_x_continuous(n.breaks=10,limits=c(min(self.sam$start)-10,max(self.sam$end)+10))+
    scale_fill_d3("category20")+guides(fill=guide_legend(title="variants"))+
    theme(text=element_text(size=15),axis.line.y = element_blank(),axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),panel.background = element_rect(fill="white"))
  if(show.entire.ref){
    p=p+scale_x_continuous(n.breaks=10)
  }
  if(return.track.info){
    return(self.sam)
  }
  if(return.snp){
    return(snp.plotdata)
  }
  return(p)
}

call.snp.from.sam<-function(sam.formated,ref){
  sam.formated$seq_without_clip=NA
  snp.plot.data=data.frame(ID=character(0),start=numeric(0),ymin=numeric(0),ymax=numeric(0),end=numeric(0),snp=character(0))
  for(i in 1:nrow(sam.formated)){
    if(grepl("S",sam.formated$CIGAR[i])){
      cigar=sam.formated$CIGAR[i]
      if(str_count(cigar,"S")==2){
        soft.clip.lengths=str_match_all(cigar,"[[:digit:]]*S")[[1]] %>% as.vector() %>% gsub("S","",.) %>% as.numeric()
        sam.formated[i,"seq_without_clip"]=substr(sam.formated$seq[i],soft.clip.lengths[1]+1,nchar(sam.formated$seq[i])-soft.clip.lengths[2])
      }else if(str_count(cigar,"S")==1&str_ends(cigar,"S")){
        soft.clip.length=str_match_all(cigar,"[[:digit:]]*S")[[1]] %>% as.vector() %>% gsub("S","",.) %>% as.numeric()
        sam.formated[i,"seq_without_clip"]=substr(sam.formated$seq[i],1,nchar(sam.formated$seq[i])-soft.clip.length)
      }else if(str_count(cigar,"S")==1&!str_ends(cigar,"S")){
        soft.clip.length=str_match_all(cigar,"[[:digit:]]*S")[[1]] %>% as.vector() %>% gsub("S","",.) %>% as.numeric()
        sam.formated[i,"seq_without_clip"]=substr(sam.formated$seq[i],soft.clip.length+1,nchar(sam.formated$seq[i]))
      }else{
        message(paste0(i,",",cigar))
        stop("Wrong CIGAR")
      }
    }else{
      sam.formated[i,"seq_without_clip"]=sam.formated$seq[i]
    }
    readseq=sam.formated$seq_without_clip[i]
    refseq=substr(ref,sam.formated$start[i],sam.formated$end[i])
    if(identical(readseq,refseq)){
      next
    }else{
      readseq=strsplit(readseq,"")[[1]]
      refseq=strsplit(refseq,"")[[1]]
      matches_for_read=sam.formated$CIGAR[i] %>% gsub("[[:digit:]]*[HS]","",.) %>% str_extract_all(.,"[[:digit:]]*[MI]")
      matches_for_read=matches_for_read[[1]]
      matches_for_ref=sam.formated$CIGAR[i] %>% gsub("[[:digit:]]*[HS]","",.) %>% str_extract_all(.,"[[:digit:]]*[MD]")
      matches_for_ref=matches_for_ref[[1]]
      if(grepl("I",sam.formated$CIGAR[i])){
        cigar_index_insertion=grep("I",matches_for_read)
        string_index_insertion=matches_for_read %>% gsub("[MI]","",.) %>% as.numeric() %>% cumsum()
        string_index_insertion=string_index_insertion[cigar_index_insertion-1]+1
        string_insertion_extend=grep("[[:digit:]]*I",matches_for_read,value=T) %>% gsub("I","",.) %>% as.numeric()
        insertion.index=mapply(FUN=seq,from=string_index_insertion,length.out=string_insertion_extend) %>% unlist()
        readseq=readseq[-insertion.index]
      }
      if(grepl("D",sam.formated$CIGAR[i])){
        cigar_index_deletion=grep("D",matches_for_ref)
        string_index_deletion=matches_for_ref %>% gsub("[MD]","",.) %>% as.numeric() %>% cumsum()
        string_index_deletion=string_index_deletion[cigar_index_deletion-1]+1
        string_deletion_extend=grep("[[:digit:]]*D",matches_for_ref,value=T) %>% gsub("D","",.) %>% as.numeric()
        deletion.index=mapply(FUN=seq,from=string_index_deletion,length.out=string_deletion_extend) %>% unlist()
        refseq[deletion.index]="-"
        readseq.mod=rep("N",length(refseq))
        readseq.mod[deletion.index]="-"
        readseq.mod[-deletion.index]=readseq
        readseq=readseq.mod
      }
#      message(i)
#      message(refseq)
#      message(readseq)
      align=compareStrings(refseq,readseq)
      mismatch.index=grep("?",align,fixed = T)
      snps=mapply(FUN=paste0,refseq[mismatch.index],rep("-",length(mismatch.index)),readseq[mismatch.index])
      snp.data=data.frame(ID=rep(sam.formated$read[i],length(mismatch.index)),start=rep(sam.formated$start[i],length(mismatch.index)),
                          ymin=rep(sam.formated$ymin[i],length(mismatch.index)),
                          ymax=rep(sam.formated$ymax[i],length(mismatch.index)))
      snp.data$start=snp.data$start+mismatch.index-1
      snp.data$end=snp.data$start
      snp.data$snp=snps
      snp.plot.data=rbind(snp.plot.data,snp.data)
    }
  }
  return(snp.plot.data)
}

SMTCoPlot<-function(srt,features,pt.size=1.1,alpha=c(0.1,6),title="nCount_SMT",log=FALSE,
                     tissue.non.tissue.colors=c("darkgrey","lightgrey"),min.count=0,max.count=NULL){
  coord=srt@images$slice1@coordinates
  coord=rownames_to_column(coord,var="barcode")
  count=FetchData(srt,vars=features) %>% as.data.frame()
  colnames(count)=colnames(count) %>% gsub("spatial_","",.)
  count=TransformRawDf(count)
  colnames(count)=c("features","barcode","count")
  if(is.null(max.count)){
    max.count=max(count$count)
  }
  count$count[count$count<min.count]=0
  count$count[count$count>max.count]=max.count
  if(log){
    count$count<-log2(count$count+1)
  }
  data=left_join(count,coord,by="barcode")
  #plot background
  plot=ggplot(data)+geom_point(data=data %>% dplyr::filter(.,tissue==0),aes(x=imagecol,y=imagerow),color=tissue.non.tissue.colors[2],size=pt.size)+
    geom_point(data=data %>% dplyr::filter(.,tissue==1),aes(x=imagecol,y=imagerow),color=tissue.non.tissue.colors[1],size=pt.size)+
    scale_y_reverse()+coord_fixed()
  #plot features
  plot=plot+geom_point(data=data %>% dplyr::filter(.,count>0),aes(x=imagecol,y=imagerow,alpha=count,color=features),size=pt.size)+
    scale_alpha_continuous(range=alpha,n.breaks=5)+scale_color_d3()
  plot=plot+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_blank(),
          plot.margin=margin(0.5,0.5,0.5,0.5,unit="cm"),
          legend.key=element_rect(fill = "white"),
          legend.text=element_text(size=12),legend.title=element_text(size=12),
          legend.margin=margin(0,0,0,0,unit="cm"),
          legend.position = "right",
          legend.direction = "vertical")+
    guides(alpha=guide_legend(title="count",override.aes=list(size=3)),
           color=guide_legend(override.aes=list(size=3)))
  return(plot)
}

TopFeatures<-function(srt,topn=10){
  data=FetchData(srt,vars=Features(srt),cells=Cells(srt))
  data=colSums(data) %>% sort(.,decreasing = T)
  return(data[1:topn])
}

corrplot.publication<-function(cors){
  corrplot(cors,method = "circle",order="hclust",type="lower",tl.col = "black",
           cl.cex = 1,tl.cex = 1,tl.srt=45,col=colorRampPalette(c("blue", "white", "red"))(100))
}

#min.count: min count of certain microbe in certain spots, the microbe's spill outs would be removed in
#all spots
filterSpillOuts<-function(srt,min.count=1000,min.overlap=5,min.R=0.5){
  DefaultAssay(srt)="Phyloseq_level"
  features=filterFeatures(srt,min.count=min.count,min.spot.occurrence=1)
  if(length(features)==0){
    message("No feature meet min.count")
    return(srt)
  }
  all.features=Features(srt)
  outs=data.frame(feat1=as.character(NULL),feat2=as.character(NULL),values=as.numeric(NULL))
  for(feature in features){
    all.features.used=all.features[!grepl(feature,all.features)]
    out=findOnevsNFeatureCor(srt=srt,feature=feature,features=all.features.used,min.coocurrence=min.overlap)
    if(nrow(out)>=1){
      outs=rbind(outs,out)
    }
  }
  outs=dplyr::filter(outs,values>=min.R)
  spillouts=outs$feat2
  valid.features=setdiff(all.features,spillouts)
  hostdata=GetAssay(srt,"Spatial")
  valid.srt=subset(srt,features=valid.features)
  valid.srt@assays$Spatial=hostdata
  return(valid.srt)
}

KL.divergence<-function(srt,feats1,feats2,norm=F){
  message("Combining data...")
  DefaultAssay(srt)="Spatial"
  feat.all=FetchData(srt,vars=Features(srt),cells=Cells(srt))
  if(any(grepl("Phyloseq_level",Assays(srt)))){
    DefaultAssay(srt)="Phyloseq_level"
    feat.microbe=FetchData(srt,vars=Features(srt),cells=Cells(srt))
    feat.all=mergeByRownames(feat.all,feat.microbe)
  }
  srt.combined.assay=CreateAssayObject(as.data.frame(t(feat.all)))
  srt.combined=srt
  srt.combined@assays$Spatial=srt.combined.assay
  message("Smoothing...")
  srt.combined=SmoothSpots(srt.combined,mode="default")
  Smoothed.data=srt.combined@assays$Smooth@counts
  pseudo.count=min(srt.combined@assays$Smooth@counts@x)*(1/5000)
  Smoothed.data=Smoothed.data+pseudo.count
  if(norm){
    Smoothed.data=t(proportions(t(Smoothed.data)))
  }
  message("filtering tissue spots...")
  tissue.spots=filterTissueSpots(srt)
  peri.spots=SelectNeighborLayers(srt=srt,spots=tissue.spots,layer_num=1)
  used.spots=union(tissue.spots,peri.spots)
  p.x=Smoothed.data[feats1,used.spots]
  q.x=Smoothed.data[feats2,used.spots]
  message("calculating matrix....")
  q.x=q.x^-1
  featin1.vs.feats2.divergences=list()
  for(feat in feats1){
    p.x.feat=t(matrix(rep(p.x[feat,],nrow(q.x)),nrow=ncol(q.x)))
    featin1.vs.feats2.divergence=p.x[feat,] %*% t(log2(p.x.feat*q.x))
    rownames(featin1.vs.feats2.divergence)=feat
    colnames(featin1.vs.feats2.divergence)=rownames(q.x)
    featin1.vs.feats2.divergences[[feat]]=featin1.vs.feats2.divergence
  }
  featin1.vs.feats2.divergences=base::Reduce(f=rbind,x=featin1.vs.feats2.divergences)
  return(featin1.vs.feats2.divergences)
}

SpatialModulePlot<-function(srt,genelist,pt.size=1.1,colors=c("#FFCC00","red"),alpha=c(0.1,6),legend.title.override=NULL,
                            min.count=NULL,max.count=NULL,legend.text.size=13,log.scale=F,legend.position="top",
                            legend.direction="horizontal",legend.pt.size=3,on.tissue.only=F){
  selfsrt=srt
#  selfsrt=try(AddModuleScore(selfsrt,list(genelist)),silent=T)
  selfsrt=try(AddModuleScore(selfsrt,genelist),silent=T)
  if(class(selfsrt)=="try-error"){
    ggplot()+geom_blank()+theme_bw()
  }else{
    SMTPlotMod(selfsrt,feature="Cluster1",pt.size=pt.size,colors=colors,alpha=alpha,
             legend.title.override=legend.title.override,min.count=min.count,max.count=max.count,
             legend.text.size=legend.text.size,log.scale=log.scale,legend.position=legend.position,
             legend.direction=legend.direction,legend.pt.size=legend.pt.size,
             on.tissue.only=on.tissue.only)
  }
}

SpatialGOPlot<-function(srt,go.id,organism,pt.size=1.1,colors=c("lightblue","red"),alpha=c(0.1,6),legend.title.override=NULL,
                        min.count=NULL,max.count=NULL,legend.text.size=13,log.scale=F,legend.position="top",
                        legend.direction="horizontal",legend.pt.size=3,on.tissue.only=T){
  if(organism=="mm"){
    organism='Mus.musculus'
  }else if(organism=="hs"){
    organism='Homo.sapians'
  }else{
    stop("invalid organism")
  }
  genes.table=try(get_anno_genes(go_ids=go.id,database=organism)$gene,silent=T)
  if(class(genes.table)=="try-error"){
    ggplot()+geom_blank()+theme_bw()
  }else{
    DefaultAssay(srt)<-"Spatial"
    SpatialModulePlot(srt=srt,genelist=list(genes.table),pt.size=pt.size,colors=colors,alpha=alpha,
                    legend.title.override=legend.title.override,min.count=min.count,max.count=max.count,
                    legend.text.size=legend.text.size,log.scale=log.scale,legend.position=legend.position,
                    legend.direction=legend.direction,legend.pt.size=legend.pt.size,
                    on.tissue.only=on.tissue.only)
  }
}

SpatialGOsubclassPlot<-function(srt,go.id,organism,pt.size=1.1,colors=c("lightblue","red"),
    alpha=c(0.1,6),min.count=NULL,max.count=NULL,legend.text.size=13,
    log.scale=F,legend.position="top",legend.direction="horizontal",legend.pt.size=3,on.tissue.only=T){
  subclass=GOfuncR::get_child_nodes(go.id) %>% dplyr::filter(.,distance==1)
  subclass.id=subclass$child_go_id
  subclass.name=subclass$child_name
  plotlist=list()
  for(i in 1:length(subclass.id)){
    message(subclass.id[i])
    plotlist[[i]]=SpatialGOPlot(srt=srt,go.id=subclass.id[i],organism=organism,pt.size=pt.size,colors=colors,alpha=alpha,
        legend.title.override=subclass.name[i],min.count=min.count,max.count=max.count,
        legend.text.size=legend.text.size,log.scale=log.scale,legend.position=legend.position,
        legend.direction=legend.direction,legend.pt.size=legend.pt.size,
        on.tissue.only=on.tissue.only)+NoLegend()+
      ggtitle(label=subclass.name[i] %>% gsub("^(.{45})(.*)$","\\1-\n\\2",.),subtitle=subclass.id[i])
  }
  ggarrange(plotlist=plotlist)
}

getGOsubclassScore<-function(srt,go.id,species="mm"){
  if(species=="mm"){
    species="Mus.musculus"
  }else if(species=="hs"){
    species="Homo.sapians"
  }
  go.ids=(get_child_nodes(go.id) %>% dplyr::filter(.,distance==1))$child_go_id
  for(id in go.ids){
    id.genes=try(get_anno_genes(id,database=species)$gene,silent = T)
    id.genes=list(id.genes)
    if(class(id.genes)=="try-error"){
      next
    }
    srt.mod=try(AddModuleScore(srt,list(id.genes),name=get_names(id)$go_name),silent=T)
    if(class(srt.mod)=="try-error"){
      next
    }
    srt=srt.mod
  }
  return(srt)
}
