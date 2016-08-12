setwd("")

## source the functions
source('haplotypeBLgenogrm.R')

## import the output file of Haploview and extract the haplotype blocks
haploblocks <- getblockhaploview(haploviewfile='example_haplo.blocks')

## after phasing your data with beagle extract the haplotypes for each block
haplos.allele <- makehaplotypes(phasedbgl='bgl.phased.gz',mapinfohap='example_haplo.map',hapblocks=haploblocks)

## to get the blocks from the above script extract them from the out list as:
## block 1...n are list 1...n
gethapalleles <- haplos.allele$HAP_ALLELES[[5]]
gethapfrequency <- haplos.allele$HAP_FREQ[[5]]

## change haplotypes within a block into
## 1: dosage/gene count == 0-HOM,1-HET,2-HOM, i.e. the sample has 'n' copies of the haplotype
haplo_dosage <- hapgenomatrix(HAP_ALLELES=haplos.allele$HAP_ALLELES[[4]],HAP_FREQ=haplos.allele$HAP_FREQ[[4]],
              MAP_info=haplos.allele$MAP_info[[4]],hapfreqThresh=0.05,outname='hap_BL1')

##### compute genomic relationship with the dosages/gene count output using vanRaden method 1 or 2
haplo_G1 <- hapGRM(haplomatrix=haplo_dosage,outputType='matrix',method ='vanRaden1',outname='hap_BL1')
haplo_G2 <- hapGRM(haplomatrix=haplo_dosage,outputType='matrix',method ='vanRaden2',outname='hap_BL1')
heatmap(haplo_G1,labRow = NA,labCol = NA)

pheno <- read.table('example_haplo.pheno',header=F)

#haplo_dosage$IID <- gsub(x=haplo_dosage$IID,pattern='^X',replacement='')
regressHAP <- merge(pheno[,-1],haplo_dosage[,-1],by=1)
regressHAP <- regressHAP[,-1]; colnames(regressHAP)[1] <- 'pheno'

summary(lm(pheno ~ ., data=regressHAP))
anova(lm(pheno ~ ., data=regressHAP))

summary(glm(pheno ~ ., family=binomial('probit'),data=regressHAP))
anova(glm(pheno ~ ., family=binomial('probit'),data=regressHAP))

for(i in 2:ncol(regressHAP)){
  print(summary(lm(pheno ~ regressHAP[,i], data=regressHAP)))
  print(summary(glm(pheno ~ regressHAP[,i], family=binomial('probit'),data=regressHAP)))
}

barplot(table(regressHAP[,c(1,2)])/sum(regressHAP$pheno),beside = T)
barplot(table(regressHAP[,c(1,3)])/sum(regressHAP$pheno),beside = T)



#system(paste('../plink.exe --chr-set 30 --nonfounders --allow-no-sex --bfile PDSalmar_chrtopsnp_bps 
#--extract snplist_top21.txt --recode 12 --out PDSalmar_haploV'))
#map <- read.table("PDSalmar_haploV.map",header=F,stringsAsFactors = F)
#write.table(map[,c(2,4)],'PDSalmar_haploV.info',quote=F,col.names=F,row.names=F)
