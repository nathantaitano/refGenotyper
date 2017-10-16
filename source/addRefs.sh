#!/bin/bash
# Uses 16 threads for alignment, and 20 gb mem. qsub command:
# qsub -q rcc-mc-30d -pe thread 16 -l mem_total=20g AddRefs.sh

tempDir="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/"
# tempDir="../temp/"
inDir="../../gbsInputData/"

outDir="../../gbsOutputData/"


refVcfHeadLine="##SAMPLE=<ID=XRef:9997,Barcode=NNNNNN,Column=12,Flowcell=NNNNNNNNNN,Flowcell_Lane=NNNNNNNNN_both,Lane=both,LibraryPrepID=9997,PlateName=A,Row=A,Sample=XRef,Status=private>"
alt1VcfHeadLine="##SAMPLE=<ID=chiltepin:9998,Barcode=NNNNNN,Column=0,Flowcell=NNNNNNNNNN,Flowcell_Lane=NNNNNNNNN_both,Lane=both,LibraryPrepID=9998,PlateName=A,Row=A,Sample=chiltepin,Status=private>"
alt2VcfHeadLine="##SAMPLE=<ID=zunla:9999,Barcode=NNNNNN,Column=0,Flowcell=NNNNNNNNNN,Flowcell_Lane=NNNNNNNNN_both,Lane=both,LibraryPrepID=9999,PlateName=A,Row=A,Sample=zunla,Status=private>"
refVcfHeadCol="XRef:9997"
alt1VcfHeadCol="chiltepin:9998"
alt2VcfHeadCol="zunla:9999"

snpContextRFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/addRefSnpContext.Rfasta"
snpContextFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/addRefSnpContext.fasta"
zunlaSaiFile=${tempDir}zunla.sai
chiltSaiFile=${tempDir}chilt.sai
zunlaSamFile=${tempDir}zunla.sam
chiltSamFile=${tempDir}chilt.sam
zunlaGenomeFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/addRef/zunla/Capsicum.annuum.L_Zunla-1_Release_2.0.fasta"
zunlaBt2Prefix="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/addRef/zunla/zunla"
chiltGenomeFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/MUMmerXRef/Capsicum.annuum.var.glabriusculum_Chiltepin_Release_2.0.fasta"
chiltBt2Prefix="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/addRef/chiltepin/chiltepin"
cm334GenomeFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/refgenome/Pepper_1.55_reform.fasta"

origVcf=${outDir}mikeyLmiss20Allhet05lmiss20imiss20.recode.vcf
vcfHeader="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/vcfHeader.txt"
rvcfFileAllGenomesAdded="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/mikeyLmiss20Allhet05lmiss20imiss20Cm334ChiltZunla.Rvcf"
vcfFileAllGenomesAdded=${outDir}mikeyLmiss20Allhet05lmiss20imiss20Cm334ChiltZunla.vcf

export PATH=${PATH}:/usr/local/R/3.3.2/bin/
export PATH=${PATH}:/usr/local/vcftools/latest/bin/


# arg 1:population pair prefix arg2:vcf arg3:number of permutations
function calcPermFst {
  p=1
  while (( $p <= $3 ))
  do
    vcftools --vcf $2 \
    --weir-fst-pop $1.${p}.pop1.txt --weir-fst-pop $1.${p}.pop2.txt \
    --out $1.perm${p}.out 2> $1.perm${p}.sum
    ((p++))
  done
}


rm $snpContextRFile # addRefsToFilteredVCF.R will append to an existing $snpContextRFile from previous runs
Rscript addRefsToFilteredVCF.R
#
fold $snpContextRFile > $snpContextFile
#

# ## align snpcontexts to other reference genomes
# ## Use Bowtie2 alignment.

############## Bowtie2 Alignment ##############
-/usr/local/bowtie2/2.2.9/bin/bowtie2 -p 16 \
  -x $chiltBt2Prefix -f \
  -U $snpContextFile \
  -S $chiltSamFile
#
/usr/local/bowtie2/2.2.9/bin/bowtie2 -p 16 \
  -x $zunlaBt2Prefix -f \
  -U $snpContextFile \
  -S $zunlaSamFile
############ End Bowtie2 Alignment ############

## parse SAM file to get SNPs for other reference genomes
Rscript addAltRefToFilteredVcf.R

## Make VCF header with new info into vcfHeader
## First, extract original vcf header
headLine=$(grep -no "#CHROM" $origVcf | cut -f 1 -d ':')
let headLess=headLine-1
head -n $headLess $origVcf > ${tempDir}origHead.vcf
echo ${refVcfHeadLine} > ${tempDir}refHead.txt
echo ${alt1VcfHeadLine} > ${tempDir}alt1Head.txt
echo ${alt2VcfHeadLine} > ${tempDir}alt2Head.txt
echo $(grep "#CHROM" $origVcf)$'\t'${refVcfHeadCol}$'\t'${alt1VcfHeadCol}$'\t'${alt2VcfHeadCol} > ${tempDir}colLine.txt
sed -i.bak "s/ /\t/g" ${tempDir}colLine.txt
cat ${tempDir}origHead.vcf ${tempDir}refHead.txt \
  ${tempDir}alt1Head.txt ${tempDir}alt2Head.txt ${tempDir}colLine.txt > $vcfHeader

cat $vcfHeader $rvcfFileAllGenomesAdded > $vcfFileAllGenomesAdded
sed -i.bak '/NA/d' $vcfFileAllGenomesAdded

vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsOutputData/fruLike_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/TusTav_keepfile.txt --weir-fst-pop ../../gbsInputData/TusTav2_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/Costeno_keepfile.txt --weir-fst-pop ../../gbsInputData/CAg_keepfile.txt \
  --out elitesFruFst

vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsOutputData/fruLike_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/TusTav_keepfile.txt --out ${outDir}fruTus1Fst
vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsOutputData/fruLike_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/TusTav2_keepfile.txt --out ${outDir}fruTus2Fst
vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsOutputData/fruLike_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/Costeno_keepfile.txt --out ${outDir}fruCostFst
vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsOutputData/fruLike_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/CAg_keepfile.txt --out ${outDir}fruCagFst

vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsInputData/TusTav_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/TusTav2_keepfile.txt --out ${outDir}tus1Tus2Fst
vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsInputData/TusTav_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/Costeno_keepfile.txt --out ${outDir}tus1CosFst
vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsInputData/TusTav_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/CAg_keepfile.txt --out ${outDir}tus1CagFst
#
vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsInputData/TusTav2_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/CAg_keepfile.txt --out tus2CagFst
vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsInputData/TusTav2_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/Costeno_keepfile.txt --out ${outDir}tus2CagFst

vcftools --vcf $vcfFileAllGenomesAdded --weir-fst-pop ../../gbsInputData/Costeno_keepfile.txt \
  --weir-fst-pop ../../gbsInputData/CAg_keepfile.txt --out ${outDir}cosCagFst

vcftools --vcf ../../gbsOutputData/mikeyLmiss20Allhet05lmiss08imiss20Cm334ChiltZunla.vcf --weir-fst-pop \
  ../../gbsOutputData/fruLike_keepfile.txt --weir-fst-pop ../../gbsInputData/Cannuum_keepfile.txt --out CfruCanFst