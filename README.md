# haplotype Analysis
Script to undertake haplotype analysis  
The scripts are R-functions developed to import haploview output and phased data (beagle v3)  
The outputs from these R script includes :  
  - extract haplotypes with a haplo-block and converting that into 1) genotype counts (0,1,2) 2) plink - ped+map file  
  - compute genomic relationship with the haplotypes (coded as genotype counts) using vanRaden method 1 or 2  

Note::  
  haplo-block = n markers make a block  
  haplotypes = the haplotypes within a block -  (maximum number of haplotypes is 2* number of samples)  

## R-scripts  
A) haplotypeBLgenogrm.R - contains 4 R-functions
    - getblockhaploview(haploviewfile)
    - makehaplotypes(phasedbgl,mapinfohap,hapblocks)
    - hapgenomatrix(HAP_ALLELES,HAP_FREQ,MAP_info,hapfreqThresh=0.05,outname)
    - hapGRM(haplomatrix,outputType,method='vanRaden1',outname)
B) runexample_haplos.R - running the example file

Example files 
A) example_haplo.ped + example_haplo.map + example_haplo.pheno - Plink --ped + --map + --pheno file   
B) example_haplo.info - haploview info file  
C) example_haplo.blocks - output of haloview (haplo-blocks)  
D) bgl.phased.gz - Beagle v3 phased data  

Procees (external)  
1. haplotypes blocks are created with haploview (using any of he block definition)  
2. haplo-blocks are exported from haploview (use the 'export options' under the 'File' tab  
3. The whole segement or a longer segment (example whole chromosome) is phased with BEAGLE v3  

Ready to use the scripts  

Data requirement  
1. The exported haplo-blocks (exported from haploview) file  
2. phased data file  
3. marker MAP information (this should be the same order and size as the marker info file used in haploview)   

See the "runexample_haplos.R" file for running the example file  

## The funtions  
- getblockhaploview(haploviewfile)  
    ... haploviewfile - the output of haloview (haplo-blocks)  

- makehaplotypes(phasedbgl,mapinfohap,hapblocks)  
    ... phasedbgl - Beagle v3 phased data  
    ... mapinfohap - marker MAP information (same order as the genotypes used in the haploview analysis)  
    ... hapblocks - the R- object of haplo-blocks extracted with the 'getblockhaploview' function above  

- hapgenomatrix(HAP_ALLELES,HAP_FREQ,MAP_info,hapfreqThresh=0.05,outname)  
    ... HAP_ALLELES - R-object that specify the haploytpes generated with the 'makehaplotypes' function above  
    ... HAP_FREQ - R-object that contains the haploytpe frequecies - also generated with 'makehaplotypes' function  
    ... hapfreqThresh - threshold for haplotypes in a block  
    ... outname - output name as plink ped+map files will be generated  

- hapGRM(haplomatrix,outputType,method='vanRaden1',outname)  
    ... haplomatrix - R-object that specify the haploytpes generated with the 'hapgenomatrix' function above  
    ... outputType - should the output be a full-matrix [use -'matrix'] or row and column wise [use-'rowcolwise']  
    ... method - vanRaden (2008) method 1 (ZZ'/sum(2pq)) [use - 'vanRaden1'] or method 2 (ZDZ'/Nsnps) [use - 'vanRaden2']  
    ... outname - output name for the grm  

