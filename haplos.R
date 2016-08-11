setwd("D:/Projects/PD_Salmar_Laura/GWAS_PRJ09/phasedata")

dat <- read.table("../PDSalmar_offspring_10pca.loco.mlma",header=T)
datchr <- dat[which(dat$Chr==21),]
thresH <- round(-log10(0.05/nrow(dat)),2)
datchr_top <- datchr[which(-log10(datchr$p)>5),]
write.table(datchr_top$SNP,'snplist_top21.txt',quote=F,col.names=F,row.names=F)

system(paste('../plink.exe --chr-set 30 --nonfounders --allow-no-sex --bfile PDSalmar_chrtopsnp_bps --extract snplist_top21.txt --recode 12 --out PDSalmar_haploV'))
map <- read.table("PDSalmar_haploV.map",header=F,stringsAsFactors = F)
write.table(map[,c(2,4)],'PDSalmar_haploV.info',quote=F,col.names=F,row.names=F)

source('getblockhaploview.R')
extractHAPs <- getblockhaploview(haploviewfile='blocks')
haplos.allele <- makehaplotypes(phasedbgl='interMS-summaryPDSalmar-phase/PDSalmar-phase.outreffile_chr21.bgl.phased.gz',
               mapinfohap='PDSalmar_haploV.map', hapblocks=extractHAPs)

gethap <- haplos.allele$HAP_ALLELES[[1]]

haplo_dosage <- hapgenomatrix(HAP_ALLELES=haplos.allele$HAP_ALLELES[[1]],HAP_FREQ=haplos.allele$HAP_FREQ[[1]],
              MAP_info=haplos.allele$MAP_info[[1]],hapfreqThresh=0.05,outname='hap_BL1')

haplo_G <- hapGRM(haplomatrix=haplo_dosage,outname='hap_BL1')


pheno <- read.table('../PDSalmar_offspring_gcta.txt',header=F)
regressHAP <- merge(pheno[,-1],haplo_dosage[,-1],by=1)
regressHAP <- regressHAP[,-1]; colnames(regressHAP)[1] <- 'pheno'

summary(lm(pheno ~ ., data=regressHAP))
summary(glm(pheno ~ ., family=binomial('probit'),data=regressHAP))
for(i in 2:ncol(regressHAP)){
  print(summary(lm(pheno ~ regressHAP[,i], data=regressHAP)))
}




barplot(table(regressHAP[,c(1,2)])/sum(regressHAP$pheno),beside = T,ylim=c(0,1))
barplot(table(regressHAP[,c(1,3)])/sum(regressHAP$pheno),beside = T,ylim=c(0,1))











X <- model.matrix(~ regressHAP[,2]-1)

for(i in 2:ncol(regressHAP)){
 a<- model.matrix(~ regressHAP[,i]-1)
if(i==2){b=a}else {b=cbind(a,b)}
}

W=data.matrix(b)

Z <- matrix(0, nrow = nrow(regressHAP), ncol = nrow(regressHAP))
diag(Z) <- 1

system(paste('../plink.exe --chr-set 30 --nonfounders --allow-no-sex --bfile ../PDSalmar_offspring --make-rel square --out PDSalmar_grm'))

G= read.table('PDSalmar_grm.rel',stringsAsFactors = F)
library(MASS)
Ginv <- ginv(data.matrix(G))
ZpZG <- crossprod(Z)+(Ginv*3.40)
LHS <- rbind(cbind(crossprod(W),crossprod(W,Z)),cbind(crossprod(Z,W),ZpZG))

y=as.matrix(data.frame(regressHAP[,1]))
RHS <- rbind(crossprod(W,y),crossprod(Z,y))

sol <- solve(LHS, RHS)

