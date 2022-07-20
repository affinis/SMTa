library('Seurat')
library('dplyr')
library('tidyverse')
library('cowplot')
library('heatmap.plus')
library("RColorBrewer")
library('ggplot2')
library("microbiomeSeq")
library("ggpubr")
library("GGally")
library("CCA")
library("ggsci")
library("phyloseq")
library("vegan")
library("ade4")
library("RVAideMemoire")
library("jpeg")
library("jsonlite")
library("igraph")
library("ggplotify")
library("parallel")
library("fpc")
library("rgl")
library('GOfuncR')

#Create Seurat Object from DBiT-seq data
CreateSrtFromDBiT<-function(mtx.path,image.path,barcode.path){
  image.dir=dirname(image.path)
  #matrix preparation
  selfmtx=read.delim(mtx.path,row.names = 1)
  selfmtx$umi=rowSums(selfmtx)
  selfmtx$B=(rownames(selfmtx) %>% strsplit(.,"x")) %>% map(.,`[[`,1) %>% unlist() %>% as.integer()
  selfmtx$A=(rownames(selfmtx) %>% strsplit(.,"x")) %>% map(.,`[[`,2) %>% unlist() %>% as.integer()
  selfbarcodelist0=read.delim(barcode.path,header=F)
  selfbarcodelist=selfbarcodelist0
  colnames(selfbarcodelist)=c("ABcombination","A","B")
  selfbarcodelist$ABcombination=selfbarcodelist$ABcombination %>% gsub("$","-1",.)
  selfmtx=left_join(selfmtx,selfbarcodelist,by=c("A","B"))
  Seuratmtx=selfmtx[,!(colnames(selfmtx) %in% c("A","B","umi"))]
  rownames(Seuratmtx)=Seuratmtx$ABcombination
  Seuratmtx$ABcombination=NULL
  Seuratmtx=Seuratmtx %>% t
  
  #craete Seurat Object
  pixelmtx_srt0=CreateSeuratObject(Seuratmtx,assay="Spatial")
  
  #create image Object for Seurat Object
  imagesize=654
  barcodelist_srt=selfbarcodelist0
  colnames(barcodelist_srt)=c("barcode","A","B")
  barcodelist_srt$barcode=barcodelist_srt$barcode %>% gsub("$","-1",.)
  barcodelist_srt$undertissue=1
  barcodelist_srt$pA=(barcodelist_srt$A-0.6)*(imagesize/15)
  barcodelist_srt$pB=(barcodelist_srt$B-0.6)*(imagesize/15)
  barcodelist_srt=barcodelist_srt[,c("barcode","undertissue","A","B","pA","pB")]
  write_csv(barcodelist_srt,file=paste0(image.dir,"/tissue_positions_list.csv"),col_names=F)
  pixelmtx_srt=pixelmtx_srt0
  system(paste0('convert -resize 654x654 ',image.path," ",paste0(image.dir,"/tissue_lowres_image.jpeg")))
  system(paste0("echo '{\"'spot_diameter_fullres'\": 7.0, \"'tissue_hires_scalef'\": 1.0, \"'fiducial_diameter_fullres'\": 25.0, \"'tissue_lowres_scalef'\": 0.3}' > ",image.dir,"/scalefactors_json.json"))
  image_srt=Read10X_JPEGImage(image.dir)
  image_srt=image_srt[Cells(pixelmtx_srt)]
  DefaultAssay(image_srt)="Spatial"
  pixelmtx_srt[["slice1"]]=image_srt
  
  return(pixelmtx_srt)
}

#perform SCT transform for srt, especially ST srt. 
SpaceWrapper<-function(sample,res){
  DefaultAssay(sample)<-"Spatial"
  sample=subset(sample,nFeature_Spatial>0)
  sample=SCTransform(sample,assay="Spatial")
  sample=RunPCA(sample, assay = "SCT", verbose = FALSE)
  sample=FindNeighbors(sample, reduction = "pca", dims = 1:30)
  sample=FindClusters(sample, verbose = FALSE,resolution=res)
  sample=RunUMAP(sample, reduction = "pca", dims = 1:30)
  return(sample)
}

TypeCluster<-function(srt,cluster,type,col.name="seurat_clusters"){
  srt@meta.data[srt@meta.data[,col.name]%in%cluster,"type"]<-type
  return(srt)
}

Selectspot<-function(data,SubsetName,...){
  chkDots(...)
  if(!exists(SubsetName,envir = globalenv())){
    assign(SubsetName,NULL,envir = globalenv())
  }
  cellist<-Cells(get(SubsetName,envir = globalenv()))
  #  cellist<-get(SubsetName,envir = globalenv())
  coords<-data@images$slice1@coordinates[,c("tissue","imagerow","imagecol")]
  p<-ggplot(coords)+geom_point(aes(imagecol,imagerow,color=tissue))+scale_y_reverse()+coord_fixed()
  cellist<-append(cellist,CellSelector(p))
  cellist<-unique(cellist)
  SeuratOb<-subset(data,cells=cellist)
  assign(SubsetName,SeuratOb,envir = globalenv())
}

#on tissue spots
filterTissueSpots<-function(tsrt,spots=NULL,counter=F){
  if(counter){
    if(is.null(spots)){
      onTissue=tsrt@images$slice1@coordinates %>% dplyr::filter(.,tissue==0) %>% rownames()
    }else{
      onTissue=spots[spots %in% (tsrt@images$slice1@coordinates %>% dplyr::filter(.,tissue==0) %>% rownames())]
    }
  }else{
    if(is.null(spots)){
      onTissue=tsrt@images$slice1@coordinates %>% dplyr::filter(.,tissue==1) %>% rownames()
    }else{
      onTissue=spots[spots %in% (tsrt@images$slice1@coordinates %>% dplyr::filter(.,tissue==1) %>% rownames())]
    }
  }
  return(onTissue)
}

#select nearest spots of certain spot
SelectNeighbors<-function(srt,spot,layout="hex"){
  ROW=srt@images$slice1@coordinates[spot,"row"]
  COL=srt@images$slice1@coordinates[spot,"col"]
  if(layout=="hex"){
    assign("upleft",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW-1,col==COL+1) %>% rownames)      
    assign("upright",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW-1,col==COL-1) %>% rownames)
    assign("downleft",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW+1,col==COL+1) %>% rownames)
    assign("downright",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW+1,col==COL-1) %>% rownames)
    assign("left",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW,col==COL-2) %>% rownames)
    assign("right",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW,col==COL+2) %>% rownames)
    return(c(upleft,upright,downleft,downright,left,right))
  }else if(layout=="cross"){
    assign("up",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW-1,col==COL) %>% rownames)
    assign("down",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW+1,col==COL) %>% rownames)
    assign("left",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW,col==COL-1) %>% rownames)
    assign("right",srt@images$slice1@coordinates %>% dplyr::filter(.,row==ROW,col==COL+1) %>% rownames)
    return(c(up,down,left,right))
  }else{
    message("Unsupported layout, 'hex' for visium, 'cross' for DBiT-Seq")
  }
}

#select nearest spots of a set of spot
SelectNeighborSets<-function(srt,spots,layout="hex"){
  PeriSpots=NULL
  for(spot in spots){
    peris=SelectNeighbors(srt,spot,layout)
    PeriSpots=append(PeriSpots,peris) %>% unique()
  }
  PeriSpots=PeriSpots[!(PeriSpots %in% spots)]
  return(PeriSpots)
}

#Select multilayers of one or a set of spots
SelectNeighborLayers<-function(srt,spots,layer_num,layout="hex"){
  selectedspots=NULL
  layer=0
  initialspots=spots
  while(layer<layer_num){
    selectedspotstmp=SelectNeighborSets(srt,spots,layout=layout)
    spots=selectedspotstmp
    selectedspots=append(selectedspots,selectedspotstmp) %>% unique()
    layer=layer+1
  }
  selectedspots=selectedspots[!(selectedspots %in% initialspots)]
  return(selectedspots)
}

#select spot to create new Seurat Object, repeat select will remove replicated spots.
CreateSeuratBySelectspot<-function(data,SubsetName,...){
  chkDots(...)
  if(!exists(SubsetName,envir=globalenv())){
    assign(SubsetName,NULL,envir = globalenv())
  }
  cellist=Cells(get(SubsetName,envir = globalenv()))
  #  cellist<-get(SubsetName,envir = globalenv())
  coords=GetTissueCoordinates(data)
  p=ggplot(coords)+geom_point(aes(imagecol,imagerow))+scale_y_reverse()
  cellist=append(cellist,CellSelector(p))
  cellist=unique(cellist)
  SeuratOb=subset(data,cells=cellist)
  assign(SubsetName,SeuratOb,envir = globalenv())
}

#a function to transform dataframe from a pure numeric matrix into a complexed one, example:
#from:
#  s1  s2  s3
#f1 0   0   1
#f2 1   2   4
#f3 2   3   2
#to:
#sample feature value
#   s1      f2    1
#   s1      f3    2
#   s2      f2    2
#   s2      f3    3
#   s3      f1    1
#   s3      f2    4
#   s3      f3    2
#note: the cutoff is the min percentage a feature should take in a sample
TransformDataframe<-function(df,cutoff){
  data=df
  cns=colnames(data)
  ret=data.frame("spot"=character(0),"features"=character(0),"values"=integer(0))
  for (i in cns){
    print(i)
    if(sum(data[,i])==0){
      next
    }
    features=(rownames(data))[data[,i]/sum(data[,i]) >= cutoff]
    if (length(features)==0){
      #      features=genus_count[,i] %>% which.max() %>% names()
      subdf=data.frame(spot=i,features="others",values=sum(data[,i]))
      ret=rbind(ret,subdf)
      next
    }
    values=data[features,i] %>% as.integer()
    spot=rep(i,length(features))
    subdf=data.frame(spot=spot,features=features,values=values)
    ret=rbind(ret,subdf)
    if (sum(data[,i])-sum(subdf$values)>0){
      otherline=data.frame(spot=i,features="others",values=sum(data[,i])-sum(subdf$values))
      ret=rbind(ret,otherline)
    }
  }
  return(ret)
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

#a function show correlation dot plot
coplot<-function(x,y,df,alpha=1,cor.coef.size=6,cor.coef.coord=c(NULL, NULL)){
  feat1=x
  feat2=y
  ggscatter(df, x = feat1, y = feat2,cor.coef.size = cor.coef.size,
            add = "reg.line", conf.int = TRUE, alpha=alpha,
            cor.coef = TRUE, cor.method = "pearson", cor.coef.coord = cor.coef.coord,
            xlab = as.character(feat1), ylab = as.character(feat2),color = "#0095DE")
}

coplot.mod<-function(df,point.color="#0095DE",cor.coef.size=3,point.size=2,
                     point.alpha=1){
  plots=list()
  void.plot=ggplot(df)+geom_blank()+
    theme(panel.background=element_rect(fill="white"),plot.margin=unit(c(0,0,0,0),"pt"))+
    xlab(NULL)+ylab(NULL)
  for(i in 1:(ncol(df)-1)){
    plot=list()
    if(i>1){
      for(m in 1:(i-1)){
        plot[[m]]=void.plot
      }
    }
    for(j in (i+1):ncol(df)){
      plot[[j-1]]=ggscatter(df,x=colnames(df)[i],y=colnames(df)[j],cor.coef=T,
          cor.coef.coord=c(min(df),max(df)*9.5/10),add="reg.line",color=point.color,
          cor.method = "pearson",cor.coef.size=cor.coef.size,size=point.size,
          conf.int = TRUE,alpha=point.alpha)+
        xlab(NULL)+ylab(NULL)+
        xlim(c(min(df),max(df)))+ylim(c(min(df),max(df)))+theme_linedraw()+
        theme(plot.margin=unit(c(0,0,0,0),"pt"),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
              text=element_text(face="bold"))
    }
    plots[[i]]=ggarrange(plotlist=plot,ncol=1,nrow=ncol(df)-1,align="v")
  }
  main=ggarrange(plotlist=plots,ncol=ncol(df)-1,nrow=1)
  xlabs.data=colnames(df)[1:(ncol(df)-1)] %>% as.data.frame()
  colnames(xlabs.data)="xlab"
  ylabs.data=colnames(df)[2:ncol(df)] %>% as.data.frame()
  colnames(ylabs.data)="ylab"
  xlabs=ggplot(xlabs.data)+geom_text(aes(x=fct_inorder(xlab),y=0.5,label=fct_inorder(xlab)))+xlab(NULL)+ylab(NULL)+
    theme(panel.background=element_rect(fill="white"),axis.text=element_blank(),axis.ticks=element_blank(),
          plot.margin=unit(c(0,0,0,0),"pt"),panel.border=element_rect(color="black",fill=NA))
  ylabs=ggplot(ylabs.data)+geom_text(aes(x=0.5,y=fct_inorder(rev(ylab)),label=fct_inorder(rev(ylab))),angle=90)+xlab(NULL)+ylab(NULL)+
    theme(panel.background=element_rect(fill="white"),axis.text=element_blank(),axis.ticks=element_blank(),
          plot.margin=unit(c(0,0,0,0),"pt"),panel.border=element_rect(color="black",fill=NA))
  ggarrange(plotlist=list(ylabs,main,void.plot,xlabs),nrow=2,ncol=2,widths=c(1,20),heights=c(20,1),align="hv")+
      theme(plot.margin=margin(0.5,0.5,0.5,0.5,"cm"))
}

#get peri-tissue spots
periTissue<-function(tsrt,layer_num,multi=T,tech="visium"){
  selfsrt=tsrt
  selflayer_num=layer_num
  selfmulti=multi
#  tissuespotidx=selfsrt@images$slice1@coordinates$tissue==1
#  tissuespot=Cells(selfsrt)[tissuespotidx]
  tissuespot=getRoughTissueEdges(selfsrt)
  if(tech=="visium"){
    selflayout="hex"
  }else if(tech=="dbit"){
    selflayout="cross"
  }else{
    message("Unsupported tech, only 'visium' and 'dbit' allowed")
  }
  res=SelectNeighborLayers(selfsrt,tissuespot,selflayer_num,layout=selflayout)
  if(selflayer_num==-1){
    outerlayer0=periTissue(selfsrt,0,tech=tech)
    res=SelectNeighborSets(selfsrt,outerlayer0,layout=selflayout)
    res=filterTissueSpots(selfsrt,res)
    res=setdiff(res,outerlayer0)
    return(res)
  }
  if(selflayer_num==0){
    outerlayer1=periTissue(selfsrt,1,tech=tech)
    res=SelectNeighborSets(selfsrt,outerlayer1,layout=selflayout)
    res=filterTissueSpots(selfsrt,res)
    return(res)
  }
  if(selfmulti==FALSE & selflayer_num>1){
    reslow=periTissue(selfsrt,selflayer_num-1,tech=tech)
    reslow=filterTissueSpots(selfsrt,reslow,T)
    res=setdiff(res,reslow)
  }
  res=filterTissueSpots(selfsrt,res,T)
  return(res)
}

extractSpotLineBinary<-function(tsrt,line_num,row=T){
  selfsrt=tsrt
  selflinenum=line_num
  if(row){
    sortedcoordbyrow=selfsrt@images$slice1@coordinates %>% arrange(.,by=row,col)
    spotlinelogic=sortedcoordbyrow$tissue[sortedcoordbyrow$row==line_num]
  }else{
    sortedcoordbycol=selfsrt@images$slice1@coordinates %>% arrange(.,by=col,row)
    spotlinelogic=sortedcoordbycol$tissue[sortedcoordbycol$col==line_num]
  }
  return(spotlinelogic)
}

extractSpotLineIdx<-function(tsrt,line_num,row=T){
  selfsrt=tsrt
  selflinenum=line_num
  if(row){
    sortedcoordbyrow=selfsrt@images$slice1@coordinates %>% arrange(.,by=row,col)
    spotlineidx=sortedcoordbyrow$col[sortedcoordbyrow$row==line_num]
  }else{
    sortedcoordbycol=selfsrt@images$slice1@coordinates %>% arrange(.,by=col,row)
    spotlineidx=sortedcoordbycol$row[sortedcoordbycol$col==line_num]
  }
  return(spotlineidx)
}

getRoughTissueEdges<-function(tsrt){
  selfsrt=tsrt
  edges=NULL
  coords=selfsrt@images$slice1@coordinates
  ROWS=coords$row %>% unique()
  COLS=coords$col %>% unique()
  for(ROW in ROWS){
    binary=extractSpotLineBinary(selfsrt,ROW)
    idx=extractSpotLineIdx(selfsrt,ROW)
    leftshiftbinary=c(binary[-1],tail(binary,1))
    rightshiftbinary=c(head(binary,1),binary[1:length(binary)-1])
    targetidx1=idx[(binary+leftshiftbinary)==1]
    targetidx2=idx[(binary+rightshiftbinary)==1]
    targetidx=c(targetidx1,targetidx2)
    barcodes=dplyr::filter(coords,row==ROW,col %in% targetidx) %>% rownames()
    edges=append(edges,barcodes)
  }
  for(COL in COLS){
    binary=extractSpotLineBinary(selfsrt,COL,F)
    idx=extractSpotLineIdx(selfsrt,COL,F)
    leftshiftbinary=c(binary[-1],tail(binary,1))
    rightshiftbinary=c(head(binary,1),binary[1:length(binary)-1])
    targetidx1=idx[(binary+leftshiftbinary)==1]
    targetidx2=idx[(binary+rightshiftbinary)==1]
    targetidx=c(targetidx1,targetidx2)
    barcodes=dplyr::filter(coords,col==COL,row %in% targetidx) %>% rownames()
    edges=append(edges,barcodes)
  }
  edges=edges %>% unique()
  edges=filterTissueSpots(selfsrt,edges)
  return(edges)
}

#SpatialDimPlot showing specified spots
viewSpots<-function(srt,spots,color="red",size=1.2){
  selfsrt=srt
  selfspots=spots
  selfcolor=c(color,"grey")
  selfsize=size
  SpatialDimPlot(selfsrt,cells.highlight=selfspots,pt.size.factor=selfsize,cols.highlight=selfcolor)
}

#create a dataframe in which first column indicate counts of lower layer while second column 
#showing counts of higher layers which are near the spots in first column, the input lowerlayer/higherlayer
#shuold be vectors composed of spot barcodes.
linkExpressionNeighbor<-function(tsrt,lowerlayer,higherlayer,metaCol="nCount_Spatial",layout="hex"){
  selfsrt=tsrt
  selflowerlayer=lowerlayer
  selfhigherlayer=higherlayer
  selfmetaCol=metaCol
  res=data.frame(layerlow=NULL,layerhigh=NULL)
  for(spot in selflowerlayer){
    neighbor=intersect(SelectNeighbors(selfsrt,spot,layout=layout),selfhigherlayer)
    nneighbor=length(neighbor)
    datalow=rep(selfsrt@meta.data[spot,selfmetaCol],nneighbor)
    datahigh=selfsrt@meta.data[neighbor,selfmetaCol]
    data=data.frame(layerlow=datalow,layerhigh=datahigh)
    res=rbind(res,data)
  }
  return(res)
}

#Diffusion plot, type can be line, bar, violin,
#layer start and end should be in c(-1,0,1,2,3,4....)
#if returndf is not null, it will use the string 
DiffusionPlot<-function(tsrt,type="violin",layerstart=-1,layerend=3,returndf=F,dfid="diffusion_df",rmOutlier=0,tech="visium"){
  selfsrt=tsrt
  selftype=type
  selflayerstart=layerstart
  selflayerend=layerend
  selfreturndf=returndf
  selfdfid=dfid
  selfdfviolin=data.frame(barcode=NULL,orig.ident=NULL,nCount_Spatial=NULL,nFeature_Spatial=NULL,layer=NULL)
  selfdfbar=data.frame(median_UMI_count=NULL,layer=NULL)
  selflineplot=list()
  for(layer in selflayerstart:selflayerend){
    layername=paste0("selflayer",layer)
    layerdfname=paste0(layername,"df")
    assign(layername,NULL)
    try(assign(layername,selfsrt@assays$Diffusion$peri.layers[[layer+2]]),silent=T)
    if(is.null(get(layername))){
      message(paste0("No corresponding info was found in Diffusion slot, calculating layer",layer," ......."))
      assign(layername,periTissue(selfsrt,layer,F,tech=tech))
    }
    assign(layerdfname,selfsrt@meta.data[get(layername),] %>% rownames_to_column(var="barcode"))
    layercolumn=rep(paste0("layer",layer),length(get(layername)))
    assign(layerdfname,cbind(get(layerdfname),layer=layercolumn))
    selfdfviolin=rbind(selfdfviolin,get(layerdfname))
    medianUMI=median((get(layerdfname))$nCount_Spatial)
    selfdfbar=rbind(selfdfbar,data.frame(median_UMI_count=medianUMI,layer=paste0("layer",layer)))
  }
  if(rmOutlier!=0){
    if(rmOutlier>0){
      i=1
      while(i<=rmOutlier){
        index=which.max(selfdfviolin$nCount_Spatial)
        selfdfviolin=selfdfviolin[-index,]
        i=i+1
      }
    }else if(rmOutlier<0){
      i=-1
      while(i>=rmOutlier){
        index=which.min(selfdfviolin$nCount_Spatial)
        selfdfviolin=selfdfviolin[-index,]
        i=i-1
      }
    }
  }
  if(selftype=="violin"){
    selfylim=max(selfdfviolin$nCount_Spatial)*1.1
    selfplot=ggplot(selfdfviolin)+
      geom_violin(aes(x=layer,y=nCount_Spatial,color=layer))+ylim(c(0,selfylim))+
      scale_color_d3(palette="category10")+
      theme_cowplot()
  }else if(selftype=="bar"){
    selfylim=max(selfdfbar$median_UMI_count)*1.1
    selfplot=ggplot(selfdfbar)+
      geom_bar(aes(x=gsub("layer","",layer),y=median_UMI_count),stat = "identity")+ylim(c(0,selfylim))+xlab("layer")+
      theme_cowplot()
  }else if(selftype=="line"){
    i=1
    for(layer in selflayerstart:(selflayerend-1)){
      layername=paste0("selflayer",layer)
      higherlayername=paste0("selflayer",layer+1)
      tmpdf=linkExpressionNeighbor(selfsrt,get(layername),get(higherlayername))
      selflineplot[[i]]=ggplot(tmpdf)+
        geom_segment(aes(x="1",xend="2",y=layerlow,yend=layerhigh))+
        ylim(c(0,max(selfsrt@meta.data$nCount_Spatial)))+xlab(NULL)+
        ylab(NULL)+theme_cowplot()+theme(axis.text.x=element_blank())
      i=i+1
      selfplot=ggarrange(plotlist=selflineplot,nrow=1,labels=paste0("layer",selflayerstart:(selflayerend-1)),vjust=-1.5)+
      theme(plot.margin=margin(1.5,0.5,0.5,0.5,"cm"))
    }
  }else if(selftype=="box"){
    selfylim=max(selfdfviolin$nCount_Spatial)*1.1
    selfplot=ggplot(selfdfviolin)+
      geom_boxplot(aes(x=layer,y=nCount_Spatial))+ylim(c(0,selfylim))+
      theme_cowplot()
  }
  if(selfreturndf){
    if(selftype=="violin" | selftype=="box"){
     selfdfviolin$sample=dfid
     return(selfdfviolin)
    }else if(selftype=="bar"){
      return(selfdfbar)
    }else{
      message("This type not supported for dataframe export..")
      stop()
    }
  }else{
    selfplot
  }
}

#load unfiltered matrix from spaceranger or cellranger
Load10X_total<-function(outdir){
  selfoutdir=outdir
  tsrt=Load10X_Spatial(selfoutdir,filter.matrix = F,filename="raw_feature_bc_matrix.h5")
  return(tsrt)
}

#get memory usage of current global environment
getGlobEnvSize<-function(){
  totalSize=0
  for (object in ls(envir=.GlobalEnv)){
    size=object.size(get(object,envir=.GlobalEnv))
    totalSize=totalSize+size
  }
  message(utils:::format.object_size(totalSize, "GB"))
}

#get peri-tissue spots by layer, and store in a new slot in Seurat object
findPeriTissue<-function(tsrt,layerstart=-1,layerend=3,tech="visium"){
  warning("This function is much slower than getPeriTissueLayers, use that function instead if layerstart>=0")
  warning("This function is going to be deprecated")
  selfsrt=tsrt
  selfstart=layerstart
  selfend=layerend
  selfsrt@assays$Diffusion=list()
  selfsrt@assays$Diffusion[["peri.layers"]]=list()
  selfsrt@assays$Diffusion[["layer.start"]]=selfstart
  selfsrt@assays$Diffusion[["layer.end"]]=selfend
  for(i in selfstart:selfend){
    message(paste0("calculating layer",i," ......"))
    index=i+2
    selfsrt@assays$Diffusion[["peri.layers"]][[index]]=periTissue(selfsrt,layer_num=i,multi=F,tech=tech)
  }
  return(selfsrt)
}


#take a look at diffusion feature of a single gene, rmOutlier is set to remove outlier data in count,
#with positive n value(like 1,2,3....) removing top n count by a increasing sorted data, with negative n,
#remove top n count by a decreasing sorted data, 0 means not remove any data.
SingleFeatureDiffusionPlot<-function(tsrt,feature,layerstart=-1,layerend=3,type="violin",rmOutlier=0,returndf=F){
  selfsrt=tsrt
  selffeature=feature
  selflayerstart=layerstart
  selflayerend=layerend
  selftype=type
  selfdfid=selffeature
  if(rmOutlier==0){
    selfOutlier=NULL
  }else{
    selfOutlier=rmOutlier
  }
  if(is.null(selfsrt@assays$Diffusion)){
    stop("Run findPeritissue() first to create diffusion slot")
  }
  if(selflayerstart<selfsrt@assays$Diffusion[["layer.start"]]
     |selflayerend>selfsrt@assays$Diffusion[["layer.end"]]){
    warning("Layer out of range, using layer data in Diffusion slot......")
    if(selflayerstart<selfsrt@assays$Diffusion[["layer.start"]]){
      selflayerstart=selfsrt@assays$Diffusion[["layer.start"]]
    }
    if(selflayerend>selfsrt@assays$Diffusion[["layer.end"]]){
      selflayerend=selfsrt@assays$Diffusion[["layer.end"]]
    }
  }
  selfbarcodes=selfsrt@assays$Diffusion$peri.layers[(selflayerstart+2):(selflayerend+2)] %>% unlist()
  selfdf=selfsrt@assays$Spatial@counts[selffeature,selfbarcodes] %>% as.data.frame() %>% rownames_to_column(.,var="barcodes")
  colnames(selfdf)<-c("barcodes",paste0(selffeature,"_nCount"))
  for(i in selflayerstart:selflayerend){
    selfdf[(selfdf$barcodes %in% selfsrt@assays$Diffusion$peri.layers[[(i+2)]]),"layer"]<-paste0("layer",i)
  }
  if(!is.null(selfOutlier)){
    if(selfOutlier>0){
      selfdf=arrange(selfdf,by=get(paste0(selffeature,"_nCount"))) %>% tail(.,(nrow(selfdf)-selfOutlier))
    }else{
      selfdf=arrange(selfdf,by=desc(get(paste0(selffeature,"_nCount")))) %>% tail(.,(nrow(selfdf)+selfOutlier))
    }
  }
  if(returndf){
    selfdf$gene=selfdfid
    return(selfdf)
  }else{
    if(selftype=="violin"){
      return(ggplot(selfdf)+
               geom_violin(aes(layer,get(paste0(selffeature,"_nCount")),color=layer))+
               theme_cowplot()+ylab(paste0(selffeature,"_nCount")))
    }else if(selftype=="box"){
      return(ggplot(selfdf)+
               geom_boxplot(aes(layer,get(paste0(selffeature,"_nCount"))))+
               theme_cowplot()+ylab(paste0(selffeature,"_nCount")))
    }
  }
}

automatedLabeling<-function(srt,featureslist=NULL){
  selfsrt=srt
  Idents(selfsrt)="seurat_clusters"
  if(is.null(featureslist)){
    selffeatures=list("B"=c("MS4A1"),
                      "Endothelial"=c("PECAM1","CD34"),
                      "Epithelial"=c("PIGR","CLDN3"),
                      "Pericyte"=c("RGS5","PDGFRB"),
                      "Fibroblast"=c("LUM","DCN","PDGFRA"),
                      "Glial"=c("S100B"),
                      "Mast"=c("TPSAB1","TPSB2","CPA3"),
                      "Myeloid"=c("C1QB","LYZ"),
                      "Plasma"=c("MZB1","IGHA1","JCHAIN"),
                      "T"=c("CD3D"),
                      "NK"=c("NCR1")
    )
  }else{
    selffeatures=featureslist
  }
  selfsrt=AddModuleScore(selfsrt,features=selffeatures)
  for(c in selfsrt@meta.data$seurat_clusters %>% sort %>% unique){
    df=(subset(selfsrt,idents=c))@meta.data[,paste0("Cluster",1:length(selffeatures))]
    max.cluster=colSums(df) %>% which.max()
    max.cluster.name=selffeatures[max.cluster] %>% names()
    selfsrt@meta.data[Cells(subset(selfsrt,idents=c)),"autotype"]=max.cluster.name
  }
  return(selfsrt)
}

getTestData<-function(sample){
  if(sample=="bc"){
    selfsrt=Load10X_Spatial(data.dir = '/data/LyuLin/Spatial/Workdirectory/10XBreastCancer/',
                           filename = 'Targeted_Visium_Human_BreastCancer_Immunology_raw_feature_bc_matrix.h5',
                           filter.matrix = F)
  }else if(sample=="crc"){
    selfsrt=Load10X_Spatial(data.dir = '/data/LyuLin/Spatial/Workdirectory/10XCRC/',
                            filename = 'Targeted_Visium_Human_ColorectalCancer_GeneSignature_raw_feature_bc_matrix.h5',
                            filter.matrix = F)
  }else if(sample=="gb"){
    selfsrt=Load10X_Spatial(data.dir = '/data/LyuLin/Spatial/Workdirectory/10XGlioblastoma/',
                            filename = 'Targeted_Visium_Human_Glioblastoma_Pan_Cancer_raw_feature_bc_matrix.h5',
                            filter.matrix = F)
  }else{
    message("Only bc(breastcancer), crc(colonrectal cancer) and gb(glioblastoma) allowed.")
  }
  return(selfsrt)
}

DiffusionDimPlot<-function(tsrt,layers=-1:3){
  selfcolors=c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B",
           "#E377C2","#7F7F7F","#BCBD22","#17BECF","#AEC7E8","#FFBB78",
           "#98DF8A","#FF9896","#C5B0D5","#C49C94","#F7B6D2","#C7C7C7")
  selfcolors=selfcolors[layers+2]
  selfcolors=append(selfcolors,"grey")
  selflayerlist=list()
  selflabels=NULL
  for(i in layers){
    if(!is.null(tsrt@assays$Diffusion$peri.layers[[(i+2)]])){
      selflayerlist[[(i+2)]]=tsrt@assays$Diffusion$peri.layers[[(i+2)]]
      selflabels=append(selflabels,paste0("layer",i))
    }else{
      message(paste0("layer ",i," data not found in Diffusion slot, will not show"))
    }
  }
  selflabels=append(selflabels,"Unselected")
#  return(selflayerlist)
#  return(selfcolors)
  SpatialDimPlot(tsrt,cells.highlight=selflayerlist,pt.size.factor=1.4)+
    scale_fill_manual(values=selfcolors,labels=selflabels)
}

getCountByLayer<-function(tsrt,layers,median=F){
  selfsrt=tsrt
  selflayers=layers
  selfdata=data.frame(barcode=NULL,nCount_Spatial=NULL,nFeature_Spatial=NULL)
  for(layer in selflayers){
    if(median){
      selfdata=DiffusionPlot(selfsrt,type="bar",layerstart=layers[1],layerend=layers[length(layers)],returndf=T)
      selfdata$layer=selfdata$layer %>% gsub("layer","",.) %>% as.numeric()
    }else{
      tmpdf=DiffusionPlot(selfsrt,type="violin",layerstart=layer,layerend=layer,returndf=T)
      tmpdf=tmpdf[,c("barcode","nCount_Spatial","nFeature_Spatial")]
      tmpdf$layer=layer
      selfdata=rbind(selfdata,tmpdf)
    }
  }
  return(selfdata)
}

batchSave<-function(pat,path){
  obs=ls(envir=.GlobalEnv) %>% grep(pat,.,value = T)
  for (object in obs){
    message(paste0("Saving ",object," ......"))
    saveRDS(object=get(object),file = paste0(path,object,".rds"))
  }
}

batchRemove<-function(pat){
  obs=ls(envir=.GlobalEnv) %>% grep(pat,.,value = T)
  for (object in obs){
    message(paste0("Deleted ",object," ......"))
    remove(object)
  }
}

mergeByRownames<-function(df1,df2){
  df<-merge(df1,df2,by="row.names",all.x=T,all.y=T)
  rownames(df)<-df$Row.names
  df$Row.names<-NULL
  return(df)
}


#a function that calculate average whitescale of each spot(one of 50x50 grid) in DBiT-Seq
#and then determine if the spot is on tissue.
getTissueStatDBiT<-function(srt,resolution=7){
  col_num=srt@images$slice1@image %>% ncol()
  row_num=srt@images$slice1@image %>% nrow()
  col_ends=(seq(0,col_num,by=col_num/50) %>% round())[-1]
  col_starts=append(1,(col_ends+1)[-length(col_ends)])
  row_ends=(seq(0,row_num,by=row_num/50) %>% round())[-1]
  row_starts=append(1,(row_ends+1)[-length(row_ends)])
  is.rgb=try(srt@images$slice1@image[,,1],silent=T)
  if('try-error' %in% class(is.rgb)){
    image.values=srt@images$slice1@image
  }else{
    image.values=srt@images$slice1@image[,,1]
  }
  for(barcode in srt@images$slice1@coordinates %>% rownames()){
    barcode_col=srt@images$slice1@coordinates[barcode,"col"]
    barcode_row=srt@images$slice1@coordinates[barcode,"row"]
    srt@images$slice1@coordinates[barcode,"grayscale"]=
    mean(image.values[row_starts[barcode_row]:row_ends[barcode_row],
                                   col_starts[barcode_col]:col_ends[barcode_col]])
  }
  kmeans=srt@images$slice1@coordinates[,"grayscale"] %>% kmeans(.,centers=resolution)
  non_tissue_cluster_num=kmeans$centers %>% which.max()
  srt@images$slice1@coordinates[,"kmeans_cluster"]=kmeans$cluster
  srt@images$slice1@coordinates[srt@images$slice1@coordinates$kmeans_cluster!=non_tissue_cluster_num,"tissue"]<-1
  srt@images$slice1@coordinates[srt@images$slice1@coordinates$kmeans_cluster==non_tissue_cluster_num,"tissue"]<-0
  return(srt)
}

Read10X_JPEGImage<-function(image.dir,filter.matrix=TRUE,...){
  image<-readJPEG(source = file.path(image.dir, "tissue_lowres_image.jpeg"))
  scale.factors<-fromJSON(txt = file.path(image.dir, "scalefactors_json.json"))
  tissue.positions<-read.csv(file = file.path(image.dir, 
              "tissue_positions_list.csv"), col.names = c("barcodes", 
              "tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, 
              as.is = TRUE, row.names = 1)
  if(filter.matrix){
    tissue.positions<-tissue.positions[which(x = tissue.positions$tissue == 
                      1), ,drop = FALSE]
  }
  unnormalized.radius<-scale.factors$fiducial_diameter_fullres* 
    scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius/max(dim(x = image))
  return(new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
        fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, 
        scale.factors$tissue_lowres_scalef), coordinates = tissue.positions, 
        spot.radius = spot.radius))
}

pixel<-function(srt,row=NULL,col=NULL,imagerow=NULL,imagecol=NULL){
  selfrow=row
  selfcol=col
  selfimagerow=imagerow
  selfimagecol=imagecol
  if(!is.null(imagerow) & !is.null(imagecol)){
    barcode=rownames(dplyr::filter(srt@images$slice1@coordinates,imagerow==selfimagerow,imagecol==selfimagecol))
  }else if(!is.null(row) & !is.null(col)){
    barcode=rownames(dplyr::filter(srt@images$slice1@coordinates,row==selfrow,col==selfcol))
  }else if(!is.null(row) & is.null(col) & is.null(imagerow) & !is.null(imagecol)){
    barcode=rownames(dplyr::filter(srt@images$slice1@coordinates,row==selfrow,imagecol==selfimagecol))
  }else if(is.null(row) & !is.null(col) & !is.null(imagerow) & is.null(imagecol)){
    barcode=rownames(dplyr::filter(srt@images$slice1@coordinates,imagerow==selfimagerow,col==selfcol))
  }else{
    message("Information not enough")
    stop()
  }
  return(barcode)
}

findImageClusters<-function(srt,clusters=5,tech="visium",smooth=F,overwrite=F){
  if(tech=="visium"){
    message("Detecting image resolution...")
    hires=srt@images$slice1@scale.factors$hires
    lowres=srt@images$slice1@scale.factors$lowres
    maxcol=srt@images$slice1@coordinates$imagecol %>% max()
    maxpixelcol=srt@images$slice1@image[,,1] %>% ncol()
    if(maxcol*hires>maxpixelcol){
      srt@images$slice1@coordinates$scaledimagerow=(srt@images$slice1@coordinates$imagerow*lowres) %>% round()
      srt@images$slice1@coordinates$scaledimagecol=(srt@images$slice1@coordinates$imagecol*lowres) %>% round()
    }else{
      srt@images$slice1@coordinates$scaledimagerow=(srt@images$slice1@coordinates$imagerow*hires) %>% round()
      srt@images$slice1@coordinates$scaledimagecol=(srt@images$slice1@coordinates$imagecol*hires) %>% round()
    }
      if(overwrite|(is.null(srt@images$slice1@coordinates$R)&is.null(srt@images$slice1@coordinates$R1))){
      message("RGB of spot center not found, determine spot center RGB...")
      for(barcode in rownames(srt@images$slice1@coordinates)){
        barcoderow=srt@images$slice1@coordinates[barcode,"scaledimagerow"]
        barcodecol=srt@images$slice1@coordinates[barcode,"scaledimagecol"]
        if(smooth){
          srt@images$slice1@coordinates[barcode,paste0("R",1:9)]=as.vector(srt@images$slice1@image[,,1][barcoderow-1:barcoderow+1,barcodecol-1:barcodecol+1])
          srt@images$slice1@coordinates[barcode,paste0("G",1:9)]=as.vector(srt@images$slice1@image[,,2][barcoderow-1:barcoderow+1,barcodecol-1:barcodecol+1])
          srt@images$slice1@coordinates[barcode,paste0("B",1:9)]=as.vector(srt@images$slice1@image[,,3][barcoderow-1:barcoderow+1,barcodecol-1:barcodecol+1])
        }else{
          srt@images$slice1@coordinates[barcode,"R"]=srt@images$slice1@image[,,1][barcoderow,barcodecol]
          srt@images$slice1@coordinates[barcode,"G"]=srt@images$slice1@image[,,2][barcoderow,barcodecol]
          srt@images$slice1@coordinates[barcode,"B"]=srt@images$slice1@image[,,3][barcoderow,barcodecol]
        }
      }
    }
    message("Finding clusters...")
    if(smooth){
      col_names=c(paste0("R",1:9),paste0("G",1:9),paste0("B",1:9))
      kMeans=kmeans(srt@images$slice1@coordinates[,col_names],centers=clusters)
    }else{
      kMeans=kmeans(srt@images$slice1@coordinates[,c("R", "G", "B")],centers=clusters)
    }
    clusters=kMeans$cluster %>% as.data.frame()
    colnames(clusters)<-"image_cluster"
    srt=AddMetaData(srt,clusters)
  }
  DefaultAssay(srt)="Spatial"
  return(srt)
}

findPeriImageCluster<-function(srt,imagecluster=NULL){
  selfsrt=srt
  selflayerstart=0
  selflayerend=1
  if(is.null(selfsrt@meta.data$image_cluster)){
    message("Image clusters not found, run findImageClusters first or define clusters by yourself.")
    stop()
  }
  if(is.null(imagecluster)){
    imagecluster=unique(selfsrt$image_cluster)
  }
  res=data.frame(barcode=NULL,parent_image_cluster=NULL,layer=NULL,image_cluster=NULL)
  for(cluster in imagecluster){
    message(paste0("Calculating image cluster ",cluster," ..."))
    selfsrt@images$slice1@coordinates[selfsrt$image_cluster==cluster,"tissue"]=1
    selfsrt@images$slice1@coordinates[selfsrt$image_cluster!=cluster,"tissue"]=0
    selfsrt=findPeriTissue(selfsrt,layerstart=selflayerstart,layerend=selflayerend,tech="visium")
    periclusterbar=selfsrt@assays$Diffusion$peri.layers %>% unlist
    periclusterbarlayers=NULL
    for(i in selflayerstart:selflayerend){
      message(paste0("layer",i))
      layerbar=selfsrt@assays$Diffusion$peri.layers[[(i+2)]]
      layer=paste0("layer",i)
      periclusterbarlayers=append(periclusterbarlayers,rep(layer,length(layerbar)))
    }
    imageclusters=selfsrt@meta.data[periclusterbar,"image_cluster"]
    clusterdf=data.frame(barcode=periclusterbar,
                         parent_image_cluster=rep(cluster,length(periclusterbar)),
                         layer=periclusterbarlayers,
                         image_cluster=imageclusters)
    res=rbind(res,clusterdf)
  }
  nonselfspots=dplyr::filter(res,layer=="layer1")
  nonselfspotsSum=plyr::count(nonselfspots,vars=c("parent_image_cluster","image_cluster"))
  imageclusterSum=selfsrt@meta.data %>% count(.,var="image_cluster")
  colnames(imageclusterSum)=c("image_cluster","total_freq")
  parentimageclusterSum=imageclusterSum
  colnames(parentimageclusterSum)=c("parent_image_cluster","parent_total_freq")
  nonselfspotsSum=left_join(nonselfspotsSum,imageclusterSum,by="image_cluster")
  nonselfspotsSum=left_join(nonselfspotsSum,parentimageclusterSum,by="parent_image_cluster")
  nonselfspotsSum$normalized.freq=(nonselfspotsSum$freq)/((nonselfspotsSum$total_freq)*(nonselfspotsSum$parent_total_freq))
  nonselfspotsSum$image_cluster=as.character(nonselfspotsSum$image_cluster)
  nonselfspotsSum$parent_image_cluster=as.character(nonselfspotsSum$parent_image_cluster)
  srt@assays$Vinicity=nonselfspotsSum
  return(srt)
}


randomPeriSpots<-function(srt,spots,max.sample.size,interval=100,repeat.times,layout="hex"){
  selflayout=layout
  totalspots=Cells(srt)
  restspots=setdiff(totalspots,spots)
  targetNeibors=SelectNeighborSets(srt,spots,layout=selflayout)
  intertargetNeibors=NULL
  for(i in seq(1,max.sample.size,interval)){
    for(repeat.count in (1:repeat.times)){
      subsample=sample(restspots,i,replace=F)
      subsampleIntargetNeibors=length(intersect(subsample,targetNeibors))
      intertargetNeibors=append(intertargetNeibors,subsampleIntargetNeibors)
    }
  }
  col.samplesize=sort(rep(seq(1,max.sample.size,interval),repeat.times))
  res=data.frame(sample.size=col.samplesize,spots.neighbor=intertargetNeibors)
  return(res)
}

SpotVinicityPlot<-function(srt){
  if(is.null(srt@assays$Vinicity)){
    message("Slot Vinicity not found in assay, run findPeriImageCluster first")
  }
  df=srt@assays$Vinicity
  links<-df[,c("parent_image_cluster","image_cluster","normalized.freq")]
  colnames(links)<-c("from","to","weight")
  links$type<-"hyperlink"
  nodes<-df[,c("parent_image_cluster","parent_image_cluster","parent_total_freq")] %>% unique()
  colnames(nodes)<-c("ID","name","freq")
  net<-graph_from_data_frame(d=links, vertices=nodes, directed=T)
  V(net)$size<-sqrt(V(net)$freq)
  V(net)$color<-pal_d3("category10")(length(V(net)$name))
  V(net)$label<-NA
  E(net)$width<-E(net)$weight*4000
  E(net)$arrow.size <-0
  edge.start<-ends(net,es=E(net), names=F)[,1]
  edge.col<-V(net)$color[edge.start]
  plot(net,edge.curved=0.2,xlim=c(-0.75,0.75),ylim=c(-0.75,0.75),edge.color=edge.col)
}

getRoughDiffusionData<-function(srt,feature.num=1000){
  inputsrt=srt
  if(is.null(inputsrt@assays$Diffusion)|is.null(inputsrt@assays$Diffusion$peri.layers[[4]])){
    message("At least from layer-1 to layer2 should be calculated, run 'findPeriTissue' first.....")
    stop()
  }
  clnum=detectCores()
  clnumUse=round(clnum*0.33)
  message(paste0(clnum," cores detected, will use ",clnumUse))
  cl=makeCluster(getOption("cl.cores", clnumUse))
  selfSingleFeatureDiffusionPlot=SingleFeatureDiffusionPlot
  getMedianCountLayer01=function(feature){
    require(dplyr)
    require(tidyverse)
    require(Seurat)
    selffeature=feature
    selfsrt=inputsrt
    rawdf=selfSingleFeatureDiffusionPlot(tsrt=selfsrt,feature=selffeature,type ="box",rmOutlier=0,returndf=T)
    mean0minus1=mean(log2(dplyr::filter(rawdf,layer %in% c("layer-1","layer0"))[,2]+1))
    mean1=mean(log2(dplyr::filter(rawdf,layer=="layer1")[,2]+1))
    mean2=mean(log2(dplyr::filter(rawdf,layer=="layer2")[,2]+1))
    df=data.frame(diffusionRatio1=mean1/mean0minus1,diffusionRatio2=mean2/mean0minus1)
    rownames(df)=selffeature
    return(df)
  }
  barcodes=inputsrt@assays$Diffusion$peri.layers[[2]]
  #"round(length(barcodes)*0.5" means the gene should occur in at least 50% layer0 spot
  ##genelist=rownames(inputsrt@assays$Spatial@counts[rowSums(inputsrt@assays$Spatial@counts[,barcodes])>=round(length(barcodes)*0.2),])
#  genelist=rownames(inputsrt)
  originset=inputsrt@assays$Spatial@counts[,barcodes] %>% rowSums() %>% sort() %>% tail(.,feature.num) %>% names()
  wholeset=inputsrt@assays$Spatial@counts %>% rowSums() %>% sort() %>% tail(.,feature.num) %>% names()
  genelist=intersect(originset,wholeset)
  message(paste0(length(genelist)," features found.."))
  message("Generating long list....")
  res=parLapply(cl,genelist,getMedianCountLayer01)
  stopCluster(cl)
  message("Merging results....")
  res=purrr::reduce(.x=res,rbind)
  res$gene=rownames(res)
  return(res)
}


RGBDimPlot<-function(srt){
  if(is.null(srt@images$slice1@coordinates$R)){
    stop("Spot center RGB not found, if this is a visium data, run findImageClusters with default args first")
  }
  data=srt@images$slice1@coordinates[,c("R","G","B")]
  data=mergeByRownames(data,srt@meta.data)
  pal=pal_d3(palette="category20")(20)
  for(i in 1:20){
    data[data$image_cluster==i,"color"]=pal[i]
  }
#  return(data)
  attach(data)
  plot3d(R,G,B,col=data$color)
  rglwidget()
}

getPeriTissueLayers<-function(srt,spots,tech="visium"){
  if(tech!="visium"){
    stop("Other tech not supported now")
  }
  clnum=detectCores()/2
  message(paste0(clnum," cores detected, will use ",clnum))
  cl=makeCluster(getOption("cl.cores", clnum))
  res=parLapply(cl,spots,CalSingleNonTissueDists,srt=srt)
  stopCluster(cl)
  res=res %>% unlist
  srt=AddMetaData(srt,res,"dist.from.tissue")
  return(srt)
}

CalSingleNonTissueDists<-function(srt,spot){
  source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')
  coords=srt@images$slice1@coordinates
  ROW1=coords %>% dplyr::filter(.,row==0)
  dist.unit=(max(ROW1$imagecol)-min(ROW1$imagecol))/63
  tissue.coords=coords %>% dplyr::filter(.,tissue==1)
  imagerow=coords[spot,"imagerow"]
  imagecol=coords[spot,"imagecol"]
  tissue.coords$row.dist=tissue.coords$imagerow-imagerow
  tissue.coords$col.dist=tissue.coords$imagecol-imagecol
  tissue.coords$dist=((tissue.coords$row.dist)^2+(tissue.coords$col.dist)^2)^(1/2)
  tissue.coords$dist.from.tissue=tissue.coords$dist/dist.unit
  layer=min(tissue.coords$dist.from.tissue)
  names(layer)=spot
  return(layer)
}

splitArea<-function(srt,split="cross"){
  coords=srt@images$slice1@coordinates
  if(split=="cross"){
    row.median=median(coords$row)
    col.median=median(coords$col)
    coords[(coords$row<row.median&coords$col<col.median),"part"]=1
    coords[(coords$row>row.median&coords$col<col.median),"part"]=3
    coords[(coords$row<row.median&coords$col>col.median),"part"]=2
    coords[(coords$row>row.median&coords$col>col.median),"part"]=4
    new.col=coords[ncol(coords)]
    srt@meta.data=mergeByRownames(srt@meta.data,new.col)
  }else{
    stop("Unsupported split method, 'cross' supported")
  }
  return(srt)
}

findSpotAggregates<-function(srt,cells,eps=2,MinPts=3){
  data=srt@images$slice1@coordinates[cells,c("row","col")]
  cluster.data=dbscan(data,eps=eps,MinPts=MinPts)
  names(cluster.data$cluster)=cells
  res=as.data.frame(cluster.data$cluster)
  colnames(res)<-"cluster"
  return(res)
}
