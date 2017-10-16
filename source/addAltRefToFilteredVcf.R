source("addRef.R")

newRefSamFile <- "/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/chilt.sam"
newRefID <- "chiltepin"
newRef2SamFile <- "/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/zunla.sam"
newRef2ID <- "zunla"
rvcfFile <- "/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/mikeyLmiss20Allhet05lmiss08imiss20Cm334.Rvcf"
newVcfFile <- "/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/mikeyLmiss20Allhet05lmiss08imiss20Cm334ChiltZunla.Rvcf"


vcf <- read.table(rvcfFile, header=T)
vcf <- addGenotypeFromSam(newRefSamFile, vcf, newRefID)
vcf <- addGenotypeFromSam(newRef2SamFile, vcf, newRef2ID)

write.table(vcf, newVcfFile, row.names = F, 
            col.names=FALSE, sep = "\t", quote = F)