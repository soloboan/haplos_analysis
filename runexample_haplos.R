setwd("")

## source the functions
source('getblockhaploview.R')

## import the output file of Haploview and extract the haplotype blocks
haploblocks <- getblockhaploview(haploviewfile='example_haplo.blocks')

## after phasing your data with beagle extract the haplotypes for each block
haplos.allele <- makehaplotypes(phasedbgl='bgl.phased.gz',mapinfohap='example_haplo.map',hapblocks=haploblocks)

## to get the blocks from the above script extract them from the out list as:
## block 1...n are list 1...n
gethapalleles <- haplos.allele$HAP_ALLELES[[1]]
gethapfrequency <- haplos.allele$HAP_FREQ[[1]]

## change haplotypes within a block into
## 1: dosage/gene count == 0-HOM,1-HET,2-HOM, i.e. the sample has 'n' copies of the haplotype
haplo_dosage <- hapgenomatrix(HAP_ALLELES=haplos.allele$HAP_ALLELES[[1]],HAP_FREQ=haplos.allele$HAP_FREQ[[1]],
              MAP_info=haplos.allele$MAP_info[[1]],hapfreqThresh=0.05,outname='hap_BL1')

##### compute genomic relationship with the dosages/gene count output using vanRaden method 1 or 2
haplo_G <- hapGRM(haplomatrix=haplo_dosage,outputType='rowcolwise',outname='hap_BL1')

pheno <- read.table('pheno.txt',header=F)
regressHAP <- merge(pheno[,-1],haplo_dosage[,-1],by=1)
regressHAP <- regressHAP[,-1]; colnames(regressHAP)[1] <- 'pheno'

summary(lm(pheno ~ ., data=regressHAP))
summary(glm(pheno ~ ., family=binomial('probit'),data=regressHAP))
for(i in 2:ncol(regressHAP)){
  print(summary(lm(pheno ~ regressHAP[,i], data=regressHAP)))
}

barplot(table(regressHAP[,c(1,2)])/sum(regressHAP$pheno),beside = T,ylim=c(0,1))
barplot(table(regressHAP[,c(1,3)])/sum(regressHAP$pheno),beside = T,ylim=c(0,1))



#system(paste('../plink.exe --chr-set 30 --nonfounders --allow-no-sex --bfile PDSalmar_chrtopsnp_bps 
#--extract snplist_top21.txt --recode 12 --out PDSalmar_haploV'))
#map <- read.table("PDSalmar_haploV.map",header=F,stringsAsFactors = F)
#write.table(map[,c(2,4)],'PDSalmar_haploV.info',quote=F,col.names=F,row.names=F)
