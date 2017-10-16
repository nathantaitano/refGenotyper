# Tools for getting the reference genome's allele, given a TASSEL-returned VCF
# (TASSEL-returned VCFs use the major allele in the REF column)
# And for creating a new VCF table that has the reference genome genotype saved as XREF
# Also writes fasta file with reference sequence around the SNP

# Features desired:
# Run with one command
# Update the reference allele (would require changing all other genotypes)
# Use per-SNP VCF format column

# getRefAllele extracts the reference allele for a certain SNP, returns 1-character genotype
# It also writes the sequence around the SNP (the two lines)
getRefAllele <- function( 
  chromPos, # a vector with 2 elements, [1] = chromosome #, [2] = pos in chromosome
  refgenomeVec, # a vector where each element is a line in the reference genome
  charsPerLine = 50, # Number of characters per line in the reference .fasta file
  chrStartLines, # a vector of line numbers in the refgenome file where the chromosome headings are located
  snpContextFile)
  {
  chrStart <- chrStartLines[chromPos[1]]
  snpLine <- chrStart + ceiling(chromPos[2] / charsPerLine)
  posInLine <- chromPos[2] %% charsPerLine
  refallele <- substr(refgenomeVec[snpLine], posInLine, posInLine)
  refContext <- paste(refgenomeVec[snpLine-1], refgenomeVec[snpLine], refgenomeVec[snpLine+1], sep="")
  posInContext <- charsPerLine + posInLine
  
  write(
    c(
      # Line 1 describes SNP chromosome/position location, 
      # and position in the reference sequence in line 2 (posHere)
      paste(c(">SNPcontext_chr",chromPos[1],"_pos",chromPos[2],
              "_posHere",posInContext), collapse = ""), 
      # Line 2 gives the reference sequence around the SNP
      refContext
    ), #SNP context data to append
    
    snpContextFile, # Filename to append data to
    
    append = T
  )
  
  return(refallele)
}

# vcf0404 <- read.table("NT_Genos_20160404_RSNPFilters_keepdepth.recode.vcf",
# comment.char = "", skip = 202, header = T)

getRefAlleles <- function(
  vcfTable, # A dataframe from VCF table, with the info lines removed (like vcf0404, above)
  refgenFile, # The path to the reference genome fasta file
  chrPrefix = ">chr",
  charsPrLine = 50, # The text that precedes each chromosome declaration in the reference fasta file
  snpContextFile = paste("addRefSnpContexts_",
                         substr(date(),9,10),substr(date(),5,7),substr(date(),21,24),
                         ".fasta",
                         sep=""))
{
  
  refgen <- readLines(refgenFile) # This step puts the whole genome into RAM. Alternatively, we could grab lines as needed using grep on linux
  chrStarts <- which(grepl(chrPrefix, refgen))
  
  refalleles <- apply(MARGIN = 1, X = vcfTable[,1:2], 
                      function(x) getRefAllele(chromPos = x, 
                                               refgenomeVec = refgen,
                                               charsPerLine = charsPrLine,
                                               chrStartLines = chrStarts,
                                               snpContextFile = snpContextFile)
                      )
  
  return(refalleles)
}

getRefEntry <- function(snp, 
                        refallele, 
                        VCFformat = c(1:5) # default is GT:AD:DP:GQ:PL. For e.g. AD:GT:DP:GQ:PL, set this as c(2,1,3,4,5)
                        ){
  GT <- c()
  ADref <- c()
  ADalt <- c()
  DP <- "5" # Not the true depth of coverage at this SNP for the reference genome. Potentially could assess this for each SNP, though it would take time. 
  GQ <- "90" # Not true GQ, like DP
  
  if(is.na(as.character(snp$REF) == as.character(refallele))){
    GTentry <- "./."
    ADalt <- "."
    ADref <- "."
    PL <- "."
    write(paste("chr",snp[1],"_pos",snp[2],"_ref",snp$REF,"_alt", collapse=""), 
          paste("refMatchesNeither_",
                substr(date(),9,10),substr(date(),5,7),substr(date(),21,24),
                ".txt",
                sep=""), # Filename to append data onto
          
          append = T
    )
  }else if(as.character(snp$REF) == as.character(refallele)){ # This is necessary because Tassel doesn't use actual ref allele (sets REF as major allele)
    GTentry <- "0/0"
    ADalt <- "0"
    ADref <- "5"
    PL <- "0,255,255" # Not true PL, like DP & GQ
  }else if(as.character(snp$ALT) == as.character(refallele)){
    GTentry <- "1/1"
    ADref <- "0"
    ADalt <- "5"
    PL <- "255,255,0"
  }else{
    GTentry <- "./."
    ADalt <- "."
    ADref <- "."
    PL <- "."
    write(paste("chr",snp[1],"_pos",snp[2],"_ref",snp$REF,"_alt", collapse=""), 
          paste("refMatchesNeither_",
                substr(date(),9,10),substr(date(),5,7),substr(date(),21,24),
                ".txt",
                sep=""), # Filename to append data onto
          
          append = T
    )
  }
  ADentry <- paste(ADref, ADalt, sep=",")
  refentry <- paste(
    c(GTentry, ADentry, DP, GQ, PL)[VCFformat], 
    collapse = ":")
  return(refentry)
}

# Requires getRefEntry (above). Could probably use apply here.
addRefsToVCF <- function(vcf, refalleles, refID="XRef"){ 
  newVCF <- cbind(vcf, rep(NA, nrow(vcf)))
  colnames(newVCF)[length(colnames(newVCF))] <- refID
  for(snp in 1:length(refalleles)){
    newVCF[snp,ncol(newVCF)] <- getRefEntry(vcf[snp,],
                                            refalleles[snp])
    
    if(snp %% 1000 == 0){ # Reports progress
      print(paste(snp, "alleles processed."))
      print(paste(snp,"th entry:", sep=""))
      print(newVCF[snp,ncol(newVCF)])
    }
  }
  return(newVCF)
}


getComp <- function(base){
  fwdBase <- c("A","C","G","T")
  compBase <- c("T","G","C","A")
  return(compBase[fwdBase==base])
}


getAlleleFromSam <- function(samRow){
  origPos <- as.numeric(str_extract(samRow[1], "(?<=posHere)[\\d]{1,3}"))
  if(as.numeric(as.character(samRow[2])) == 16 & samRow[6] == "150M"){
    allele <- getComp(substr(samRow[10],150-(origPos-1),150-(origPos-1)))
  }else if(as.numeric(as.character(samRow[2])) == 0 & samRow[6] == "150M"){
    allele <- substr(as.character(samRow[10]), origPos, origPos)
  }else{ 
    allele <- "N"
  }
  return(allele)
}

addGenotypeFromSam <- function(samFile, vcfIn, refID){
  library(stringr)
  
  samDf <- read.table(samFile, comment.char = "@", fill = T) # Cols are sometimes inconsistent past the first 10, but we only use 1,6,10 anyway.
  samDf[,1] <- as.character(samDf[,1])
  samDf[,10] <- as.character(samDf[,10])
  samDf[,6] <- as.character(samDf[,6])
  
  samAlleles <- apply(MARGIN = 1, X = samDf, FUN = getAlleleFromSam)
  
  snpIDs <- paste("S", str_extract(samDf[,1],"(?<=chr)\\d{1,2}"),
                  "_", str_extract(samDf[,1],"(?<=pos)\\d+"),
                  sep="")
  snpCHR <- as.numeric(str_extract(samDf[,1],"(?<=chr)\\d{1,2}"))
  snpPOS <- as.numeric(str_extract(samDf[,1],"(?<=pos)\\d+"))
  samAlleles <- samAlleles[order(snpCHR,snpPOS)]
  
  vcfDf <- addRefsToVCF(vcfIn, samAlleles, refID=refID)
  
  for(i in 100:ncol(vcfDf)){
    vcfDf[,i] <- as.character(vcfDf[,i])
  }
  
  return(vcfDf)
}

addFromSnpCallGenome <- function(genomeFile, vcf){
  genomeAlleles <- getRefAlleles(vcf, genomeFile)
  newVcf <- addRefsToVCF(vcf, genomeAlleles)
  return(newVcf)
}

getRvcfFromVcf <- function(origVcfFile){
  first100kLines <- readLines(origVcfFile, n=100000)
  startLine <- which(str_detect(first100kLines, "#CHROM"))
  vcf <- read.table(origVcfFile,
                    comment.char = "", skip = startLine-1, header = T)
  return(vcf)
}

#############################################################################################
####### MUMmer XRef Section (for show-aligns output) ########################################
#############################################################################################

# Takes output from the show-aligns program in MUMmer
# Returns a table giving information, including error count, start, and end, for each alignment
makeAlnTable <- function(rawAlns){
  beginLineNos <- which(grepl("BEGIN", rawAlns))
  endLineNos <- which(grepl("END", rawAlns))
  
  refAlignStartPoses <- as.numeric(str_extract(rawAlns[beginLineNos],
                                               "(?<=\\[ [+-][1] )\\d+"))
  refAlignEndPoses <- as.numeric(str_extract(rawAlns[beginLineNos],
                                             "(?<=\\d - )\\d+(?= | )"))
  querAlignStartPoses <- as.numeric(str_extract(rawAlns[beginLineNos],
                                                "(?<=\\| [+-][1] )\\d+"))
  querAlignEndPoses <- as.numeric(str_extract(rawAlns[beginLineNos],
                                              "(?<=\\d - )\\d+(?= \\])"))
  
  alnLen <- refAlignEndPoses - refAlignStartPoses
  
  numErrs <- c()
  for(i in 1:length(beginLineNos)){
    ithAln <- rawAlns[beginLineNos[i]:endLineNos[i]]
    ithAnn <- ithAln[seq(6, length(ithAln), 4)]
    justAnn <- sapply(ithAnn, function(x) substr(x, 12, 49+12))
    justAnn <- paste(justAnn, collapse = "")
    numErrs <- c(numErrs, 
                 length(str_extract_all(justAnn,"\\^+")[[1]])
    )
  }
  alnTable <- data.frame(refStartPos = refAlignStartPoses, querStartPos = querAlignStartPoses,
                         refEndPos = refAlignEndPoses, querEndPos = querAlignEndPoses,
                         rawBeginLine = beginLineNos, rawEndLine = endLineNos,
                         errCnt = numErrs, alnLen = alnLen)
  return(alnTable)
}

#####The two functions below really could have been combined into one more efficient step#######
# Takes the table produced above, a reference position (numeric), and show-aligns output
# Returns two lines containing the best alignment surrounding the reference position provided
getQuerTargLine <- function(refPos, rawAlns, alnTable){
  matchRows <- which(alnTable$refStartPos <= refPos & alnTable$refEndPos >= refPos)
  matchErrRates <- alnTable[matchRows,"errCnt"] / alnTable[matchRows, "alnLen"]
  if(length(matchErrRates)==0){
    return(NA)
  }
  matchRow <- matchRows[which(matchErrRates == min(matchErrRates))]
  if(length(matchRow) > 1){
    matchRow <- sample(matchRow, 1) # picks a random match in case of error rate ties
  }
  
  rawMatch <- rawAlns[alnTable[matchRow, "rawBeginLine"]:alnTable[matchRow, "rawEndLine"]]
  sideDex <- as.numeric(str_extract(rawMatch, "^\\d+(?=\\s)"))
  sideDex <- sideDex[-which(is.na(sideDex))]
  refSideDex <- sideDex[seq(1,length(sideDex),2)]
  refLineDex <- unlist(sapply(refSideDex, 
                              function(x) which(str_detect(
                                rawMatch, paste(as.character(x), 
                                                "\\s+[cagt\\.]{1,49}$", sep="")))))
  querSideDex <- sideDex[seq(2,length(sideDex),2)]
  querLineDex <- unlist(sapply(querSideDex, 
                               function(x) which(str_detect(
                                 rawMatch, paste(as.character(x), 
                                                 "\\s+[cagt\\.]{1,49}$", sep="")))))
  
  for(i in 1:length(refSideDex)){
    if(i < length(refSideDex)){
      if(refSideDex[i] <= refPos & refSideDex[i+1] > refPos){
        refTargLine <- rawMatch[refLineDex[i]]
        querTargLine <- rawMatch[refLineDex[i]+1]
      }
    }else if(i == length(refSideDex) & refPos >= refSideDex[i]){
      refTargLine <- rawMatch[refLineDex[i]]
      querTargLine <- rawMatch[refLineDex[i]+1]
    }
  }
  
  return(list(refTargLine = refTargLine, querTargLine = querTargLine))
}

# Takes the single-line alignment selected from show-aligns output, and the reference position
# Returns the query genome allele at the alignment to that reference position
getQuerAllele <- function(refPos, querTargLine, refTargLine){
  targDex <- as.numeric(str_extract(refTargLine, "\\d+"))
  targBases <- strsplit(str_extract(refTargLine, "[acgtn\\.]+"), "")[[1]]
  querTargBases <- strsplit(str_extract(querTargLine, "[acgtn\\.]+"), "")[[1]]
  targPoses <- c()
  for(i in 1:length(targBases)){
    if(targBases[i] == "."){
      targPoses <- c(targPoses, NA)
    }else{
      targPoses <- c(targPoses, targDex + i-1 - sum(is.na(targPoses)))
    }
  }
  querBase <- querTargBases[which(targPoses == refPos)]
  querBaseOptions <- c("a","c","g","t","n",".")
  outBaseOptions <- c("A","C","G","T","N","N")
  return(outBaseOptions[querBaseOptions == querBase])
}

#Takes one chromosome of snp positions at a time (as numeric vector), original vcf, refID and...
# the output of MUMmer's show-aligns function
#Returns a vcf with the aligned query genome's genotypes added to the vcf file
addMumAlignedToVCF <- function(rawAlnFile, refPoses, vcf, refID){
  print("Loading alignment file...")
  rawAlns <- readLines(rawAlnFile)
  print("Alignment file loaded. Making table of alignments...")
  alnTable <- makeAlnTable(rawAlns)
  print("Alignment table complete. Extracting alleles...")
  querAlleles <- c()
  for(pos in refPoses){
    targLines <- getQuerTargLine(pos, rawAlns, alnTable)
    
    if(is.na(targLines)){
      querAlleles <- c(querAlleles, "N")
      print(paste("Allele", querAlleles[length(querAlleles)], "extracted from position:",pos))
      next
    }
    
    querAlleles <- c(querAlleles, 
                     getQuerAllele(refPos = pos, 
                                   querTargLine = targLines$querTargLine,
                                   refTargLine = targLines$refTargLine))
    if(length(querAlleles) %% 500 == 0){
      print(paste("Allele", querAlleles[length(querAlleles)], "extracted from position:",pos))
    }
  }
  
  newVCF <- addRefsToVCF(vcf=vcf, refalleles=querAlleles, refID = refID)
  return(newVCF)
}