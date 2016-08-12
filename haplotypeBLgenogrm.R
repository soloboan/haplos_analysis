getblockhaploview <- function(haploviewfile){
  hapblocks <- read.table(haploviewfile,fill=T,stringsAsFactors=F,na.strings=" ")
  hapblocks <- hapblocks[which(hapblocks$V1=='BLOCK'),]
  hapblocks <- cbind.data.frame(NR=1:nrow(hapblocks),hapblocks)
  extractHAPs <- data.frame(blocknr=character(),hapstr=numeric(),hapend=numeric())
  for(i in 1:nrow(hapblocks)){
    gethapblock <- hapblocks[i,5:ncol(hapblocks)]
    gethapblock <- na.omit(t(gethapblock))
    gethapblock <- as.vector(gethapblock[which(gethapblock[,1]!=""),])
    hapblockstr <- as.numeric(as.vector(gsub("!","",gethapblock[1])))
    hapblockend <- as.numeric(as.vector(gsub("!","",gethapblock[length(gethapblock)])))
    hapsinfo <- cbind.data.frame(blocknr=hapblocks[i,3],hapstr=hapblockstr,hapend=hapblockend)
    extractHAPs <- rbind.data.frame(extractHAPs,hapsinfo)
  }
  return(extractHAPs)
}

makehaplotypes<- function(phasedbgl,mapinfohap,hapblocks){
  hapgeno.bgl <- read.table(paste(phasedbgl),header=T,stringsAsFactors = F)
  map <- read.table(paste(mapinfohap),header=F,stringsAsFactors = F)
  snpids <- hapgeno.bgl[,2]
  hapgeno.bgl <- t(hapgeno.bgl[,-1:-2])
  colnames(hapgeno.bgl) <- snpids
  stringhaps <-function(x){paste(x,sep="",collapse="")}
  nanim <- nrow(hapgeno.bgl)/2
  nhaps <- nrow(hapgeno.bgl)
  animids <- rownames(hapgeno.bgl)[seq(1,nrow(hapgeno.bgl),2)]
  
  for(j in 1:nrow(hapblocks)){
    snpids <- map[c(hapblocks[j,2]:hapblocks[j,3]),2]
    ### mat haplo
    dathaplo <- hapgeno.bgl[seq(1,nrow(hapgeno.bgl),2),snpids]
    mat.haps <- data.frame(mathaplo=apply(dathaplo,1,stringhaps))
    
    ### pat haplo
    dathaplo <- hapgeno.bgl[seq(2,nrow(hapgeno.bgl),2),snpids]
    pat.haps <- data.frame(pathaplo=apply(dathaplo,1,stringhaps))
    
    # haplotypes
    haps <- cbind.data.frame(pat.haps,mat.haps)
    numhaps <- unclass(c(as.factor(haps[,1]),as.factor(haps[,2])))
    numhaps <- unclass(as.factor(c(as.vector(haps[,1]),as.vector(haps[,2]))))
    
    haps$pathaps.num <- numhaps[1:nanim]
    haps$mathaps.num <- numhaps[(1+nanim):nhaps]
    haps$FID <- animids
    haps$IID <- animids
    haps <- haps[,c(5:6,3:4,1:2)]
    
    namehaplo <- paste("hapBlock-",j,sep="")
    phase.haplonum <- rbind.data.frame(data.frame(HAP=haps$pathaps.num),data.frame(HAP=haps$mathaps.num))
    phase.haplo <- rbind.data.frame(data.frame(HAP=haps$pathaplo),data.frame(HAP=haps$mathaplo))
    hapnumalleles <- unique(cbind.data.frame(HAPNUM=phase.haplonum,HAPALLELE=phase.haplo))
    hapnumalleles <- hapnumalleles[order(hapnumalleles$HAP),]
    
    HAPs.results <- data.frame(table(phase.haplonum))
    HAPs.results$propHAP <- round(HAPs.results$Freq/sum(HAPs.results$Freq),4)
    HAPs.results <- merge(HAPs.results,hapnumalleles,by=1)
    HAPs.results <- HAPs.results[order(-HAPs.results$Freq),]
    
    maphaplo <- map[map[,2] %in% snpids,]
    
    haps <- list(haps); names(haps) <- namehaplo
    HAPs.results <- list(HAPs.results); names(HAPs.results) <- namehaplo
    maphaplo <- list(maphaplo); names(maphaplo) <- namehaplo
    
    if(j==1){
      haplos <- haps
      haplos.freq <- HAPs.results
      maphaplotype <- maphaplo
    } else {
      haplos<-c(haplos,haps)
      haplos.freq<-c(haplos.freq,HAPs.results)
      maphaplotype<-c(maphaplotype,maphaplo)
    }
  }
  return(list(HAP_ALLELES=haplos,HAP_FREQ=haplos.freq,MAP_info=maphaplotype))
}


hapgenomatrix <- function(HAP_ALLELES,HAP_FREQ,MAP_info,hapfreqThresh=0.05,outname){
  freqthaps <- HAP_FREQ[which(HAP_FREQ$propHAP>=hapfreqThresh),]
  hapclasses <- as.vector(freqthaps$phase.haplonum)
  nhap <- as.numeric(length(hapclasses))
  nanim=nrow(HAP_ALLELES)
  
  haplomatrix <- matrix(0,nrow=nanim,ncol=nhap)
  names_HAP <- paste('HAP_',hapclasses,sep='')
  colnames(haplomatrix) <- names_HAP
  for(j in 1:nrow(HAP_ALLELES)){
    for (i in 1:length(hapclasses)){
      pat=ifelse(HAP_ALLELES[j,'pathaps.num']==hapclasses[i],1,0)
      mat=ifelse(HAP_ALLELES[j,'mathaps.num']==hapclasses[i],1,0)
      haplomatrix[j,i] <- pat+mat
    }
  }
  haplomatrix <- cbind.data.frame(HAP_ALLELES[,1:2],haplomatrix)
  
  convert2ped <- function(dat){
    dat <- data.frame(dat)
    dat[which(dat[,1]==0),1] <- '11'; dat[which(dat[,1]==1),1] <- '12'; dat[which(dat[,1]==2),1] <- '22'
    return(dat)
  }
  haps2ped <- data.frame(apply(haplomatrix[,-1:-2],2,convert2ped))
  colnames(haps2ped) <- names_HAP
  haps2ped <- cbind.data.frame(HAP_ALLELES[,1:2],haps2ped)
  
  write.table(haps2ped,paste(outname,'.ped',sep=''),quote=F,col.names=F,row.names=F)
  mapR <- cbind.data.frame(CHR=unique(MAP_info[,1]),name=names_HAP,GP=0,PP=1:length(names_HAP))
  write.table(mapR,paste(outname,'.map',sep=''),quote=F,col.names=F,row.names=F)
  
  return(haplomatrix)
}

hapGRM <- function(haplomatrix,method='vanRaden1',outputType='rowcolwise',outname){
  animids <- as.vector(haplomatrix[,2])
  M=data.matrix(haplomatrix[,-1:-2])
  if(method=='vanRaden1'){
    Z=scale(M,center=T,scale=F)
    ZZp <- tcrossprod(Z)
    p=colMeans(M)/2
    K=2*sum(p*(1-p))
    G=ZZp/K
    #hist(diag(G),breaks=100)
    #hist(G[lower.tri(G,diag=F)],breaks=100)
  } else if (method=='vanRaden2'){
    Z=scale(M,center=T,scale=T)
    ZZp <- tcrossprod(Z)
    K=ncol(M)
    G=ZZp/K
    #hist(diag(G),breaks=100)
    #hist(G[lower.tri(G,diag=F)],breaks=100)
  }
  colnames(G) <- animids
  rownames(G) <- animids
  if(outputType=='rowcolwise'){
  require(reshape2)
  rowiseG <- melt(G)
  idsnum <- data.frame(IID=unclass(as.factor(c(rowiseG[,1],rowiseG[,2]))))
  idsnumL <- nrow(idsnum)
  idsnumbind <- cbind(idsnum[1:(idsnumL/2),],idsnum[(1+(idsnumL/2)):idsnumL,])
  rowiseG <- cbind(idsnumbind,rowiseG)
  write.table(rowiseG,paste(outname,'.grm',sep=''),quote=F,col.names=F,row.names=F)
  return(rowiseG)
  } else {
  write.table(G,paste(outname,'.grm',sep=''),quote=F,col.names=F,row.names=F)
  return(G)
  }
}
