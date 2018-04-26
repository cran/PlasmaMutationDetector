# release 1.3   : first release fully running
# release 1.4   : in BuildCtrlErrorRate, hotspot file was not used correctly -> corrected
#                 in DetectPlasmaMutation, adjbox is now computed on the finite values of the 250 largest values
# release 1.4.1 : correct dependencies for R version 3.3.1 and CRAN submission
# release 1.4.2 : add data for examples
# release 1.4.3 : reduce data size and change PrepareLibrary function to not use Exac file
# release 1.4.4 : reduce data size more to fit CRAN policy (< 5MB) and add possibility to use extension .bam.bai for bai files in MAF_from_BAM
# release 1.4.5 : call to data(...) corrected with data(...,envir = environment())
# release 1.4.6 : LICENSE correction
# release 1.4.7 : remove "don't run" with shrinked examples
# release 1.4.8 : remove warnings coming from search of InDel with less than 3 bases + from geom_point and geom_text (see YVES 5/09/2016 in MAF_from_BAM)
# release 1.4.9 : add file.path() to example to run under windows
# release 1.5.0 : remove file.path() and check whether pathnames end properly by a '/' (see # TO RUN UNDER WINDOWS)
# release 1.5.1 : downsample example to fit the <5s time constraint of the CRAN
# release 1.5.2 : use Nicolas example to fit the <5s time constraint of the CRAN and add the possibility to use zip BER
# release 1.6   : adapt the code to accept no large deletion, remove one mistake and simplify code
# release 1.6.1 : bug introduced during simplifcation fixed
# release 1.6.2 : bug introduced during simplifcation fixed
# release 1.6.3 : in BuildCtrlErrorRate, use hotspot file as data.frame 27/11/2017. Now hotspot file should contains a variable name chrpos on the first line
# release 1.6.4 : save file in tempdir() by default to satisfy CRAN policy e.g. for debian distribution
# release 1.6.5 : debug 1.6.4
# release 1.7.0 : take care of mutations showing allele frequency larger than 50% // clean the code

# #' @importFrom GenomicRanges GRanges
#' @importClassesFrom S4Vectors FilterRules
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment rowRanges
#' @import GenomicRanges
#' @import VariantAnnotation
#' @importFrom Rsamtools PileupParam BamFile scanBam pileup filterBam ScanBamParam scanBamWhat
#' @import ggplot2
#' @import grid
#' @importFrom graphics abline plot points text
#' @importFrom stats dnorm integrate median pbinom quantile binom.test
#' @importFrom utils data download.file read.delim read.table write.table unzip
#'
#' @encoding{utf-8}

# Bioconductor packages
# Install :
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
# biocLite("VariantAnnotation")
# biocLite('Rsamtools')
# biocLite('rtracklayer')
# biocLite("SummarizedExperiment")

# R packages
require(robustbase)
require(ggplot2)
require(grid)
require(GenomicRanges)
require(VariantAnnotation)
require(Rsamtools)
require(rtracklayer)

##################### UTILITY FUNCTIONS #####################

c2n = function(x) as.numeric(as.character(x))

sbscore = function(refplus, refmin, altplus, altmin) {

  if (is.na(altplus)|is.na(altmin)) return(NA)

  rmM = function(x,y) min(x,y)/max(x,y)
  z = refplus*altmin/(refmin*altplus)
  log( rmM(refplus,refmin)/rmM(altplus,altmin) * (z+1/z) )
}

minus.logit = function(P) {
  sapply(P,function(x) {if (!is.nan(x)&(x<1)) return(-log(x/(1-x))) else return(-Inf)})
  }

pvalue = function(E,N=NULL,E0=NULL,N0=NULL,with.info=FALSE,use.int=TRUE) {
  # return the pvalue
  if (is.null(N)) {N=E[2]; E0=E[3]; N0=E[4]; E=E[1]}
  if (E0==0) E0=1/sqrt(N0) # provides an upper bound for the p-value by assuming control error rate equals 1/N0^(3/2)
  p0.hat = E0/N0
  sigma0.hat = sqrt(p0.hat*(1-p0.hat)/N0)
  if (with.info) {
    print(paste0('p0.hat=',p0.hat,' --- p.hat=',E/N))
    print(sigma0.hat)
  }
  psi = function(y) pbinom(E,N,pmax(pmin(p0.hat-sigma0.hat*y,1),0),lower.tail=F)*dnorm(y)
  if (use.int) p.val = integrate(psi,lower=-30,upper=30,subdivisions=10000L)$value else p.val = binom.test(E,N,p0.hat,'great')$p.value
  p.val
}


find.outliers = function(P,chr.pos,thr.R=6.5,thr.P=0.000001,with.info=FALSE) {

  P[which(P==1)]=0.9999999999 ### just to have finite values
  P[which(P==0)]=0.0000000001 ### just to have finite values
  nP = length(P)
  names(P) = chr.pos
  L.tri = sort(minus.logit(P))
  dL = diff(L.tri);
  L.right = L.tri[-1];
  L.left = L.tri[-nP];
  sL.sqrt=sign(L.left)*(abs(L.left))^0.5 # to take the square root without changing the sign
  Ratio = dL/sL.sqrt
  outliers = which( Ratio>thr.R & L.right>minus.logit(thr.P))

  ### ADD YVES 29/03/2016
  # if (length(outliers)>0) outliers = which(L.right>=min(L.right[min(outliers)]))

  ind = names(outliers)

  if (with.info) { ### just for debugging purpose. Not used.
    plot(L.tri[-1][L.left>0],Ratio[L.left>0],col=1)
    abline(v=minus.logit(thr.P),col='red'); abline(h=thr.R,col='green'); abline(h=6,col='cyan'); abline(h=7,col='violet');
    # if (length(outliers)>0) outliers = which(L.right>=min(L.right[min(outliers)]))
    if (length(outliers)>0) points(L.right[outliers],Ratio[outliers],col='red',pch='x')
    if (length(outliers)>0) text(L.right[outliers],Ratio[outliers],labels=ind,cex=0.6,pos=2,col=1)
    if (length(outliers)>0) print(sort(P,decreasing=T)[outliers+1])
  }

  ind
}


pileupFreq = function(pileupres) {
  nucleotides = levels(pileupres$nucleotide)
  res = split(pileupres, pileupres$seqnames)
  res = lapply(res, function (x) split(x, x$pos))
  res = lapply(res, function (positionsplit) {
    nuctab = lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tablecounts = sapply(nucleotides, function (n)
        sum(each$count[each$nucleotide == n]))
      tablecountsplus = sapply(nucleotides, function (n)
        sum(each$count[each$nucleotide == n & each$strand=="+"]))
      tablecountsminus = sapply(nucleotides, function (n)
        sum(each$count[each$nucleotide == n & each$strand=="-"]))
      strandbias = sapply(nucleotides, function (n)
        max(0,each$count[each$strand=="+" & each$nucleotide == n])/
          (max(0,each$count[each$strand=="+" & each$nucleotide == n]) + max(0,each$count[each$strand=="-" & each$nucleotide == n]))
      )
      c(chr,pos, tablecounts, tablecountsplus, tablecountsminus, strandbias)
    })
    nuctab = data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) = NULL
    nuctab
  })
  res = data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) = NULL
  colnames(res) = c("seqnames","start",levels(pileupres$nucleotide),
                    paste0(levels(pileupres$nucleotide),"_plus"),
                    paste0(levels(pileupres$nucleotide),"_minus"),
                    paste0(levels(pileupres$nucleotide),"_SB"))
  res[3:ncol(res)] = apply(res[3:ncol(res)], 2, as.numeric)
  res
}

############ END OF UTILITY FUNCTIONS ##################


############ DOCUMENTED FUNCTIONS ######################

#' function PrepareLibrary
#'
#' Define the Genomic Ranges and Genomic Positions covered by the AmpliSeq™ Panel to include in the study and define SNP positions to exclude from the study.
#' Trimming amplicon ends is performed if specified. This function is mostly useful if you want to add some SNP positions which are not existing in the
#' positions_ranges.rda file provided within the package. It is provided to be able to reconstruct \code{positions_ranges.rda} data.
#'
#' @param info.dir, char, name of the folder containing the library information files (default 'Info/')
#' @param bed.filename, char, name of a BED table (tab-delimited) describing the Panel (with first 3 columns: "chr" (ex:chr1), "start position" (ex:115252190), "end position" (ex:115252305), i.e. the Ion AmpliSeq™ Colon and Lung Cancer Research Panel v2 (default 'lungcolonV2.bed.txt' as provided in the inst/extdata/Info folder of the package).
#' @param snp.filename, char, name of the vcf file describing known SNP positions, obtained from ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz (default 'ExAC.r0.3.sites.vep.vcf.gz'). It requires a corresponding TBI file to be in the same folder (obtained from ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz.tbi)
#' @param snp.extra, a vector of char, a vector of extra known snp positions manually curated (ex:"chrN:XXXXXXXXX")
#' @param output.name, char, filename to save \code{pos_ind} and \code{pos_snp} (default 'positions_ranges.rda')
#' @param output.dir, char, directory where to save \code{pos_ind} and \code{pos_snp} (default \code{info.dir})
#' @param load.from.broad.insitute, boolean, if TRUE load \code{snp.filename} from Broad Institute ftp server otherwise use the file positions_ranges_broad.rda (default FALSE)
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @seealso positions_ranges,
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
#'
#' @return Save the following variables in a .rda file defined by \code{output.name} in the folder defined by \code{output.dir}:
#' \itemize{
#' \item \code{pos_ranges}, a GRanges descriptor of amplicon positions
#' \item \code{pos_ind}, a vector of char "chrN:XXXXXXXXX", defining ALL index positions
#' \item \code{pos_snp}, a vector of char "chrN:XXXXXXXXX", defining SNP positions
#' }
#'
#' @export PrepareLibrary
#'
#' @examples
#'    bad.pos = "chr7:15478"
#'    PrepareLibrary(info.dir='./',snp.extra=bad.pos,output.dir=paste0(tempdir(),'/'))
#'
#'
PrepareLibrary = function(info.dir='Info/',
                          bed.filename='lungcolonV2.bed.txt', # amplicon positions
                          snp.filename='ExAC.r0.3.sites.vep.vcf.gz', # known snp
                          snp.extra=c("chr4:1807909","chr7:140481511", "chr18:48586344", "chr14:105246474", "chr19:1223055"), # extra snp
                          output.name='positions_ranges.rda',
                          output.dir=info.dir,
                          load.from.broad.insitute=FALSE
) {

  IRanges = NULL

  # Find SNP
  snp.file = paste0(info.dir,snp.filename)
  tbi.file = paste0(snp.file,'.tbi')
  if (!file.exists(snp.file) | !file.exists(tbi.file)) {
    print(paste0('File ',snp.file,' or ',tbi.file,' do not exist in ',info.dir,' ...'))
    if (load.from.broad.insitute) {
      print('Download ExAC release 0.3 from Broad Institute ... please wait !!!')
      if (!file.exists(snp.file)) download.file(paste0('ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/',snp.filename), snp.file, method='libcurl')
      if (!file.exists(tbi.file)) download.file(paste0('ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/',snp.filename,'.tbi'), paste0(snp.file,'.tbi'), method='libcurl')
    } else {
      # DEFAULT : use provided file positions_ranges.rda
      print('Use positions_ranges.rda of the package')
      data(positions_ranges,envir = environment())
    }
  }

  if (file.exists(snp.file) & file.exists(tbi.file)) {
    # Read bed file
    bed.file = paste0(info.dir,bed.filename)
    if (!file.exists(bed.file)) stop(paste0(bed.file,' not found in ',info.dir))

    bed = read.delim(bed.file, header=FALSE) # original amplicons
    # Build index of all positions
    pos_ind = NULL
    for(i in 1:nrow(bed)) pos_ind = c(pos_ind,paste0(bed[i,1],':',bed[i,2]:bed[i,3]))

    # original positions for each amplicons
    pos_ranges = GRanges(bed[,1], IRanges(bed[,2], bed[,3]))

    # Build our snp selection
    vcf = readVcf(snp.file, "hg19", param=ScanVcfParam(which=GRanges(gsub('chr','',bed[,1]), IRanges(bed[,2], bed[,3]))))
    snp = cbind(info(vcf),as.data.frame(rowRanges(vcf))) # extract vcf information
    snp$AF = sapply(snp$AF, max) # consider only max allele frequency
    snp = snp[which(snp$AF>1e-4),] # only "real" snp with Allele Frequency >1e-4
    pos_snp = paste0('chr',snp$seqnames,':',snp$start) # build snp positions
  }

  pos_snp = c(pos_snp,snp.extra) # add the SNP positions you want


  output.name = paste0(output.dir,output.name)
  save(pos_ranges,pos_ind,pos_snp,file=output.name)
  print('File saved in ...')
  return(output.name)
}


#' function MAF_from_BAM
#'
#' Read BAM files  and create MAF file. BAMfiles are stored in a sub-folder '/rBAM'.
#' MAF files are intermediate files stored in a sub-folder '/BER'.
#' MAF files contain the raw counts of A,T,C,G, insertion, deletion, insertion>2bp, deletion >2bp for strand plus and stand minus.
#' Note : we strongly recommand to externally recalibrate BAM files using tools like GATK.
#'
#' @param study.dir, char, name of the folder containing the rBAM directory  (default 'Plasma/'). The typical folder hierarchy will consist of 'Plasma/rBAM'
#' @param input.filenames, a vector of char (default NULL), the names of the BAM files to process. If NULL all BAM files in the rBAM folder will be processed
#' @param bai.ext, char, filename extension of the bai files (default '.bai')
#' @param pos_ranges.file, char, name of the Rdata file containing the three variables \code{pos_ind}, \code{pos_snp} and \code{pos_ranges} as build by the function \code{PrepareLibrary}. Default NULL, use the position_ranges.rda provided, used for our analysis.
#' @param force, boolean, (default FALSE) if TRUE force all computations to all files including already processed ones
#' @param output.dir, char, name of the folder to save results  (default \code{study.dir}).
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
#'
#' @return the path/names of the MAF files
#'
#' @export MAF_from_BAM
#'
#' @examples
#'   \dontrun{
#'      ctrl.dir = system.file("extdata", "4test_only/ctrl/", package = "PlasmaMutationDetector")
#'      if (substr(ctrl.dir,nchar(ctrl.dir),nchar(ctrl.dir))!='/')
#'        ctrl.dir = paste0(ctrl.dir,'/') # TO RUN UNDER WINDOWS
#'      MAF_from_BAM(ctrl.dir,force=TRUE,output.dir=paste0(tempdir(),'/'))
#'    }
#'
MAF_from_BAM = function(study.dir='Plasma/',input.filenames=NULL,bai.ext='.bai',
                        pos_ranges.file=NULL,force=FALSE,output.dir=study.dir) {

  pos_ind = pos_snp = pos_ranges = NULL

  if (is.null(pos_ranges.file)) {
    data(positions_ranges,envir=environment())
  } else {
    if (!file.exists(pos_ranges.file)) stop(paste('STOP :::',pos_ranges.file,' does not exist'))
    load(pos_ranges.file)
  }

  ######## Input/Output folders
  bam.dir = paste0(study.dir,'rBAM/')        # Folder for BAM files and BAI index files
  ber.dir = paste0(study.dir,'BER/')         # Default folder for BER files
  ber.outdir = paste0(output.dir,'BER/')     # Folder to save Background Error Rate files extracted from BAM files
  if (!dir.exists(ber.outdir)) dir.create(ber.outdir)

  ######## pileup parameters
  pileup.param = PileupParam(max_depth=100000, min_base_quality=20, min_mapq=5, min_nucleotide_depth=1, min_minor_allele_depth=0, # min_base_quality=20, min_mapq=5
                             distinguish_strands=TRUE, distinguish_nucleotides=TRUE, ignore_query_Ns=TRUE, include_deletions=TRUE,
                             include_insertions=TRUE)

  ######## Samples to process
  if (is.null(input.filenames)) {
    print(paste0("Looking for bam and bai in ",bam.dir,' ...'))
    input.filenames = dir(bam.dir, pattern="*.bam$") # list all BAM files in bam.dir
  }

  bam.files = paste0(bam.dir,input.filenames) # build bam.files as paste0(study.dir,'rBAM/',input.filenames[i])

  if (length(bam.files)==0) stop(paste('STOP ::: no bam file found ... NOTE: bam files need to be in a subdirectory called rBAM'))

  ber.files = rep(NA,length(bam.files))

  for( g in 1:length(bam.files)){

    bam.name = bam.files[g]
    bai.name = gsub('.bam$',bai.ext,bam.name)
    ber.name = gsub(bam.dir,ber.dir,gsub('.bam$','_MAF.txt',bam.name),fixed=TRUE) # local BER
    berzip.name = paste0(ber.name,'.zip') # local zip BER
    ber.outname = gsub(bam.dir,ber.outdir,gsub('.bam$','_MAF.txt',bam.name),fixed=TRUE) # BER to create
    berzip.outname = paste0(ber.outname,'.zip') # BER to create in zip version

    ber.files[g] = ber.outname

    if (!force & file.exists(ber.name)) {
      ber.files[g] = ber.name
      next()
    }
    if (!force & file.exists(ber.outname)) {
      ber.files[g] = ber.outname
      next()
    }

    if (!force & file.exists(berzip.name)) {
      ber.files[g] = unzip(zipfile=berzip.name,exdir=ber.outdir)
      # write.table(readLines(unz(berzip.name,rev(strsplit(ber.name,'/')[[1]])[1])),file=ber.outname,sep=' ',row.names=FALSE,col.names=FALSE,quote=FALSE)
      next()
    }
    if (!force & file.exists(berzip.outname)) {
      ber.files[g] = ber.files[g] = unzip(zipfile=berzip.outname,exdir=ber.outdir)
      # write.table(readLines(unz(berzip.outname,rev(strsplit(ber.name,'/')[[1]])[1])),file=ber.outname,sep=' ',row.names=FALSE,col.names=FALSE,quote=FALSE)
      next()
    }

    print(paste0('Processing MAF from file ',bam.name))

    ### Load BAM and BAI files
    if (!file.exists(bam.name)) stop(paste('STOP :::',bam.name,'not found ... NOTE : bam files need to be in a subdirectory called rBAM'))
    if (!file.exists(bai.name)) {
      bai.sndname = paste0(bam.name,'.bai')
      if (!file.exists(bai.sndname)) stop(paste('STOP :::',bai.name,' and ',bai.sndname,'not found ... NOTE : bai files need to be in a subdirectory called rBAM'))
      bai.name = bai.sndname
    }
    bf = BamFile(bam.name, index=bai.name)

    ### Import, count, index, filter, sort, and merge BAM
    ft = scanBam(bf, scanBamParam=ScanBamParam(which=pos_ranges, what=scanBamWhat()), pileupParam=pileup.param)

    ### Count SNV from the BAM file
    res = pileup(bf, scanBamParam=ScanBamParam(which=pos_ranges, what=scanBamWhat()), pileupParam=pileup.param)

    ### Call function pileupFreq to get a data.frame of counts
    freq = pileupFreq(res)

    freq$chrpos = paste(freq[,1],freq[,2], sep=":") # positions existing in the BAM file
    select = freq[which(freq$chrpos %in% pos_ind),] # extract positions of interest


    ### Compute coverture; major coverture, major frequency, background frequency, second major coverture
    select[,'cov']=select[,3]+select[,4]+select[,5]+select[,6]
    select[,'maj']=apply(select[,3:6],1,max)
    select[,'fmaj']=select$maj/select$cov
    select[,'bruitdefond']=1-select$fmaj
    select[,'MA'] = apply(select[,3:6],1,sort)[3,]
    # select[,'M3'] = apply(select[,3:6],1,sort)[2,]
    # select[,'M4'] = apply(select[,3:6],1,sort)[1,]

    ### Select and count reads with DELETIONS >2 nucleotides
    dels = grep("[D]", ft[[1]][["cigar"]], invert=FALSE)  # index reads carrying deletions
    numdel = regmatches(as.character(ft[[1]]$cigar[dels]), gregexpr('([0-9]+)D',ft[[1]]$cigar[dels])) # get number of deleted bp
    numdel.max = sapply(lapply(numdel,function(x) as.numeric(gsub('D','',x)) ),max) # for reads with multiple deletions, get the maximum number of deleted bp
    selection = dels[which(as.numeric(numdel.max)>2)] # index reads carrying deletions > 2 bp
    subft = sapply(ft[[1]],function(x) x[selection] ) # generate filter for temporary BAMfile export that contain reads with deletion >2bp only

    ### Select read ids with deletion >2bp
    filter = S4Vectors::FilterRules(list(KeepQname = function(ft) ft$qname %in% subft$qname))
    ### Build a temporary BAM file with selected reads
    dest = filterBam(bam.name, tempfile(),index=bai.name, filter=filter, param=ScanBamParam(what="qname"))
    ### Import this temporary BAM file
    gal = rtracklayer::import(dest,format='bam', param=ScanBamParam(what=(c("qname", "seq", "qual"))) ) # import the temporary BAMfile that contain reads with deletion >2bp only

    ### Select read ids with long deletions, skip the deletions < 3bp that could co-occur in the same read by replacing CIGAR D by CIGAR N
    for (z in 1:length(gal@cigar)) {
      mm = regmatches(gal@cigar[z], gregexpr('[0-9]+[A-Z]',gal@cigar[z]))
      ind = grep('D',mm[[1]])       # YVES 5/09/2016
      if (length(ind)>0) {          # YVES 5/09/2016
        jnd = which(as.numeric(gsub('D','',mm[[1]][ind]))<3)                    # YVES 5/09/2016
        if (length(jnd)>0) mm[[1]][ind[jnd]] = gsub('D','N',mm[[1]][ind[jnd]])  # YVES 5/09/2016
      }                             # YVES 5/09/2016
      # mm[[1]][which(as.numeric(as.character(gsub("D", "", mm[[1]]))) < 3)]=gsub("D", "N", mm[[1]][which(as.numeric(as.character(gsub("D", "", mm[[1]]))) < 3)])
      gal@cigar[z] = paste(mm[[1]], collapse="")
    }
    ### Build a temporary BAM file with selected reads
    dest = rtracklayer::export(gal, BamFile(tempfile())) # export the temporary BAMfile cleared of any deletions < 3bp
    ### Import this temporary BAM file
    res = pileup(dest, scanBamParam=ScanBamParam(which=pos_ranges, what=scanBamWhat()), pileupParam=pileup.param)

    ### Build a data.frame of long deletion counts
    if (nrow(res)>0) {
      freq = pileupFreq(res)	#read the cleared temporary BAMfile
      contig = freq[which(freq[,'-'] > 0),c(1:9,17,25,33)] #create dataframe for long deletions
      contig$contig = rep(NA,nrow(contig)) #add a new variable to define starting point of deletion
      # print(names(contig))

      # YVES 08/11/2017
      # A simpler way to find large deletions (>2bp) and assign : YVES 8/11/2017
      # frequency consistence changes in this version
      loc = dloc = 1
      while(loc<nrow(contig)) {
        while ((loc+dloc<=nrow(contig)) & (as.numeric(contig$start[loc])+dloc==as.numeric(contig$start[loc+dloc])) &
               (abs(mean(contig[loc+(0:(dloc-1)),'-'])-contig[loc+dloc,'-'])/mean(contig[loc+(0:(dloc-1)),'-']) < 0.2)) dloc = dloc+1
        dloc = dloc-1 # now dloc is the number of bases with same deletion as loc
        if (dloc>1) { # found large deletion (>2bp)
          # print(loc+(0:dloc))
          # cat(contig$seqnames[loc+(0:dloc)],'\t'); cat(contig$seqnames[loc+dloc+(1:3)],'\n')
          # cat(contig$start[loc+(0:dloc)],'\t'); cat(contig$start[loc+dloc+(1:3)],'\n')
          # cat((contig[loc+(1:dloc),'-']-contig[loc,'-'])/contig[loc,'-'],'\t'); cat((contig[loc+dloc+(1:3),'-']-contig[loc,'-'])/contig[loc,'-'],'\n')
          contig$contig[loc+(0:dloc)] = loc
          loc = loc+dloc+1
        } else {
          loc = loc+1
        }

        dloc = 1
      }

      ### remove positions with inconsistent deletions (<2bp or different frequency)
      contig2 = contig[!is.na(contig$contig),]

      contig2$chrpos = rep(NA,nrow(contig2)) # to avoid pb when empty

      if (nrow(contig2)>0) contig2$chrpos=paste0(contig2$seqnames,":",contig2$start)

      # YVES 08/11/2017
      ### A warning to prevent the false duplication of deletions
      if (any(duplicated(contig2$chrpos)))
          cat(' !!!!!!! warning  !!!!! duplicates in deletion for sample',bam.files[g],'\n')

    } else {
      contig2 = as.data.frame(matrix(nrow=0,ncol=14))
      colnames(contig2)=c("seqnames","start","A","C","G","T","N","=","-","-_plus","-_minus","-_SB","contig","chrpos")
    }

    ### Select and count reads with INSERTIONS >2 nucleotides
    ins = grep("[I]", ft[[1]][["cigar"]], invert=FALSE) # index reads carrying insertions
    numins = regmatches(as.character(ft[[1]]$cigar[ins]), gregexpr('([0-9]+)I',ft[[1]]$cigar[ins])) # get number of inserted bp
    numins2 = sapply(sapply(numins,function(x){gsub('I','',x)}),max) # for reads with multiple insertions, get the maximum number of inserted bp
    selection=ins[which(as.numeric(numins2)>2)] # index reads carrying deletions > 2bp
    subft = sapply(ft[[1]],function(x) x[selection])

    ### Select read ids with deletion >2bp only
    filter = S4Vectors::FilterRules(list(KeepQname = function(ft) ft$qname %in% subft$qname))
    ### Build a temporary BAM file with selected reads
    dest = filterBam(bam.name, tempfile(),index=bai.name, filter=filter, param=ScanBamParam(what="qname"))
    ### Import this temporary BAM file
    gal = rtracklayer::import(dest,format='bam', param=ScanBamParam(what=(c("qname", "seq", "qual"))) ) # import the temporary BAMfile that contain reads with deletion >2bp only

    # For reads with long insertions, skip the insertions < 3bp that could co-occur in the same read by replacing CIGAR I by CIGAR S
    for (z in 1:length(gal@cigar)){
      mm = regmatches(gal@cigar[z], gregexpr('[0-9]+[A-Z]',gal@cigar[z]))
      ind = grep('I',mm[[1]])       # YVES 5/09/2016
      if (length(ind)>0) {          # YVES 5/09/2016
        jnd = which(as.numeric(gsub('I','',mm[[1]][ind]))<3)                    # YVES 5/09/2016
        if (length(jnd)>0) mm[[1]][ind[jnd]] = gsub('I','S',mm[[1]][ind[jnd]])  # YVES 5/09/2016
      }                             # YVES 5/09/2016
      # mm[[1]][which(as.numeric(as.character(gsub("I", "", mm[[1]]))) < 3)]=gsub("I", "S", mm[[1]][which(as.numeric(as.character(gsub("I", "", mm[[1]]))) < 3)])
      gal@cigar[z] = paste(mm[[1]], collapse="")
    }
    ### Build a temporary BAM file with selected reads
    dest = rtracklayer::export(gal, BamFile(tempfile())) # export the temporary BAMfile cleared of any insertsions < 3bp
    ### Import this temporary BAM file
    res = pileup(dest, scanBamParam=ScanBamParam(which=pos_ranges, what=scanBamWhat()), pileupParam=pileup.param)

    if(nrow(res)>0) {
      freq = pileupFreq(res)		# read the cleared temporary BAMfile
      incontig = freq[which(freq[,'+'] > 0),c(1:8,10,18,26,34)] # create dataframe for long insertions
    } else {
      incontig = as.data.frame(matrix(nrow=0,ncol=13))
      colnames(incontig)=c("seqnames","start","A","C","G","T","N","=","+","+_plus","+_minus","+_SB","chrpos")
    }
    unlink(dest)

    ##### Merge dataframes for SNV, DELETIONS>2bp and INSERTIONS>2bp
    if (nrow(incontig)>0) incontig[,'chrpos']=paste(incontig[,'seqnames'],":",incontig[,'start'],sep='')

    mergeindels = merge(contig2[,c(9:14)], incontig[,9:13], by='chrpos', all=TRUE)
    colnames(mergeindels) = c('chrpos','longdel','deletion_plus','deletion_minus','longdel_SB','contig','longins','insertion_plus','insertion_minus','longins_SB')
    select = merge(select, mergeindels, by='chrpos', all=TRUE)
    select[is.na(select$longdel),'longdel'] = 0
    select[is.na(select$longins),'longins'] = 0
    select[,'longINDEL'] = apply(select[,c('longdel','longins')],1,max)

    ### save ber.file
    write.table(select[!duplicated(select$chrpos),], file=ber.outname, row.names=FALSE, col.names=TRUE) ####  save temporary files
  }

  return(ber.files)

}


#' function BuildCtrlErrorRate
#'
#' Compute the SNV Position-Error Rates and INDEL Position-Error Rates from control samples (available in the control directory \code{ctrl.dir}).
#' This function requires MAF files, that will be automatically generated if not present in the specified control folder.
#' SNV PER is computed as the sum in control samples of SNV background counts / sum in control samples of depths where SNV background counts = depth - major allele count.
#' INDEL PER is computed as sum in control samples of INDEL background counts / sum in control samples of depths where INDEL background counts = sum of insertion and deletion counts.
#'
#' @param ctrl.dir, char, foldername containing the control files (default 'Plasma ctrl/'). The typical folder hierarchy will consist of 'Plasma ctrl/rBAM'
#' @param bai.ext, char, filename extension of the bai files (default '.bai')
#' @param pos_ranges.file, char, name of the Rdata file containing the three variables \code{pos_ind}, \code{pos_snp} and \code{pos_ranges} as build by the function \code{PrepareLibrary}. Default NULL, use the position_ranges.rda provided, used for our analysis.
#' @param hotspot.file, char, name of the text file containing a list of the genomic positions of the hotspots (default NULL, read the provide hotspot.txt, see \code{hotspot})
#' @param force, boolean, (default FALSE) if TRUE force all computations to all files including already processed ones
#' @param output.dir, char, name of the folder to save results (default \code{ctrl.dir}).
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
#'
#' @return the number of processed files
#'
#' @export BuildCtrlErrorRate
#'
#' @examples
#' \dontrun{
#'    ctrl.dir = system.file("extdata", "4test_only/ctrl/", package = "PlasmaMutationDetector")
#'    if (substr(ctrl.dir,nchar(ctrl.dir),nchar(ctrl.dir))!='/')
#'      ctrl.dir = paste0(ctrl.dir,'/') # TO RUN UNDER WINDOWS
#'    BuildCtrlErrorRate(ctrl.dir,output.dir=paste0(tempdir(),'/'))
#'    }
#'
BuildCtrlErrorRate = function(ctrl.dir='Plasma ctrl/',bai.ext='.bai',pos_ranges.file=NULL,hotspot.file=NULL,
                              force=FALSE,output.dir=ctrl.dir) {


  ### read position informations
  pos_snp = pos_ind = NULL
  if (is.null(pos_ranges.file)) {
    data(positions_ranges,envir=environment())
  } else {
    if (!file.exists(pos_ranges.file)) stop(paste('STOP :::',pos_ranges.file,'does not exist'))
    load(pos_ranges.file)
  }

  ### read hotspot positions
  hotspot = NULL
  if (is.null(hotspot.file)) {
    data(hotspot,envir=environment())
  } else {
    if (!file.exists(hotspot.file)) stop(paste('STOP :::',hotspot.file,'does not exist'))
    hotspot = read.delim(hotspot.file, header=TRUE) # YVES 27/11/2017 FALSE -> TRUE
    if (!any(colnames(hotspot)=='chrpos')) stop(paste('STOP :::',hotspot.file,'has to contain a variable named chrpos (defined on first line)'))
  }

  ### extract counts from all BAM files found in ctrl.dir
  ctrl.files =  MAF_from_BAM(study.dir=ctrl.dir,input.filenames=NULL,bai.ext=bai.ext,pos_ranges.file=pos_ranges.file,force=force,output.dir=output.dir)

  print(ctrl.files)

  # ### list Background-Error-Rate files
  # ctrl.files = dir(paste0(ctrl.dir,'BER'), pattern="*.txt$", full.names=TRUE) # list all Background Error Rate files # before 1.6.5

  ### Assemble background errors from control samples
  X = read.table(ctrl.files[1], header=TRUE, sep=" ")

  ### utility function
  ATGC.ind = which(colnames(X) %in% c('A','T','G','C'))
  ATGC = colnames(X)[ATGC.ind]
  fct.MAJ = function(X) {
    data.frame(chrpos=as.character(X$chrpos),REF=apply(X,1,function(x) ATGC[which(x[ATGC.ind]==x['maj'])[1]]))
  }

  # extract main variant = genome ref
  MAJ = fct.MAJ(X)

  if (length(ctrl.files)>1) ind.g = 2:length(ctrl.files) else ind.g = 1 # just for testing purpose with one bam
  for (g in ind.g) {
    Y = read.table(ctrl.files[g],header=T,sep=" ")
    X = merge(Y,X,by="chrpos",all=TRUE,suffixes=c(as.character(g),''))
    MAJ = merge(fct.MAJ(Y),MAJ,by='chrpos',all=TRUE,suffixes=c(as.character(g),''))
  }

  ok = (X$chrpos %in% pos_ind)
  not.ok = (X$chrpos %in% pos_snp)
  X = X[ok,] #### ensure to select positions in defined genomic ranges
  X = X[!not.ok,] #### ensure to remove SNP positions and background error rate for SNPs (E0) and INDELS (E0indel) in samples
  MAJ = MAJ[ok,]
  MAJ = MAJ[!not.ok,]
  # print(X[X$chrpos=='chr6:170530538',])
  # print(X[X$chrpos=='chr10:39080921',])

  X.col = colnames(X)
  col.cov = grep('cov',X.col) # column indexes of the coverture
  col.del = grep('longdel($|[^_])',X.col) # column indexes of the longdel
  col.ins = grep('longins($|[^_])',X.col) # column indexes of the longins
  col.maj = grep('^maj',X.col) # column indexes of the major

  #### Check for good coverage
  # not.good.coverage = apply(X[,col.cov],1,function(x) which(is.na(x)))
  # print(not.good.coverage[1:3])
  # for (i in 1:length(not.good.coverage))
  #   if (length(not.good.coverage[[i]])>0) {
  #     print(paste('PROBLEM at positions',X$chrpos[i],'in MAF/BER :'))
  #     print(not.good.coverage[[i]])
  #     print(ctrl.files[not.good.coverage[[i]]])
  #   }

  sum.no.na = function(x) sum(x,na.rm=TRUE)
  q95 = function(x) quantile(x,0.95,na.rm=TRUE)

  N0 = apply(X[,col.cov],1,sum.no.na) # Sum of the depths of all control samples for each position
  E0 = apply(X[,col.cov]-X[,col.maj],1,sum.no.na) # Sum of the SNV background counts of all control samples for each position

  X = X[which(!is.na(N0) | N0 > 5000),]  # Remove genomic positions which are not covered in all control samples or with sum of depths < 5,000x
  MAJ = MAJ[which(!is.na(N0) | N0 > 5000),]

  REF = apply(MAJ,1,function(x) names(which.max(table(as.character(x[-1])))))
  MAJ$the.ref = REF
  write.table(MAJ, file='MAJ.txt', row.names=FALSE, col.names=TRUE)

  ### SNV rates
  p.sain = E0/N0 # Control SNV Position-Error Rate
  up.sain = apply((X[,col.cov]-X[,col.maj])/X[,col.cov],1,q95) # 95th percentile for Control SNV Position-Error Rate

  ### INDEL rates
  E0indel = apply(cbind(X[,col.del], X[,col.ins]),1,sum.no.na) # Sum of the INDEL background counts of all control samples for each position
  indel.p.sain = E0indel/N0 # Control INDEL Position-Error Rate
  indel.up.sain = apply((X[,col.del] + X[,col.ins]) / X[,col.cov], 1, q95) # 95th percentile for Control INDEL Position-Error Rate

  # Joint dataframe containing Control SNV and INDEL Position-Error Rate together with their 95th percentile, and the raw counts data
  chrpos = as.character(X$chrpos)
  background_error_rate = data.frame(cbind(chrpos,N0,E0,p.sain,up.sain,E0indel,indel.p.sain,indel.up.sain,REF))

  # Precise what are the hotspot positions
  background_error_rate$hotspot = 'Non-hotspot'
  background_error_rate[which(background_error_rate$chrpos %in% hotspot$chrpos ),'hotspot'] = 'Hotspot'

  # Save Background error file...
  backround_error_file = paste0(output.dir,"background_error_rate.txt")
  print('Backround error rate saved in ...')
  print(paste0('... ',backround_error_file))
  write.table(background_error_rate, file=backround_error_file, row.names=FALSE, col.names=TRUE)

  return(length(ctrl.files))
}


#' function LoadBackgroundErrorRate
#'
#' This function will load the background error rates created from the controls using the function BuildCtrlErrorRate
#'
#' @param pos_ranges.file, char, name of the Rdata file containing the three variables \code{pos_ind}, \code{pos_snp}, \code{pos_ranges} as build by the function \code{PrepareLibrary}. Default NULL, use the position_ranges.rda provides that we used for our analysis.
#' @param ber.ctrl.file, char, pathname of the file providing the background error rates obtained from the controls (default NULL use the provided background error rates obtained from our 29 controls). See \code{background_error_rate.txt} data and \code{BuildCtrlErrorRate} function.
#' @param n.trim, integer, number of base positions trimmed at the ends of each amplicon (default 8)
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
#'
#' @return the adapted background error rate
#'
LoadBackgroundErrorRate = function(pos_ranges.file,ber.ctrl.file,n.trim) {
  ### read background error rate
  background_error_rate = NULL
  if (is.null(ber.ctrl.file)) {
    data('background_error_rate',envir=environment()) # load default background error rate obtained from the 29 controls
  } else {
    if (!file.exists(ber.ctrl.file)) stop(paste('STOP:::File',ber.ctrl.file,'does not exist'))
    background_error_rate = read.table(ber.ctrl.file, header=T, sep=" ")
  }

  ### read position informations
  pos_ranges = pos_snp = pos_ind = NULL
  if (is.null(pos_ranges.file)) {
    data(positions_ranges,envir=environment())
  } else {
    if (!file.exists(pos_ranges.file)) stop(paste('STOP :::',pos_ranges.file,' does not exist'))
    load(pos_ranges.file)
  }

  if (n.trim>0) {
    pos_ind = NULL
    chrpos = rep(pos_ranges@seqnames@values,times=pos_ranges@seqnames@lengths)
    for(i in 1:length(chrpos)) {
      start = pos_ranges@ranges@start[i]
      end = start+pos_ranges@ranges@width[i]-1
      pos_ind = c(pos_ind,paste0(chrpos[i],':',(start+n.trim):(end-n.trim)))
    }
  }

  background_error_rate = subset(background_error_rate,chrpos %in% pos_ind)
  background_error_rate
}

#' function UpperPart
#'
#' Find the locations inside a set of values higher than an upper quantile of these values
#'
#' @param values, numeric, the set of values
#' @param q, the upper quantile
#'
#' @return the indices of the highest values, larger than the given quantile
#'
UpperPart = function(values,q) {
  noisy = which(values > quantile(values,q))
  noisy
}


#' function DetectPlasmaMutation
#'
#' This is the main function of the package that calls mutations by comparing at each genomic position the SNV or INDEL frequencies computed in one tested sample to
#' the SNV or INDEL Position-Error Rates computed from several control samples by a binomial test. An outlier detection is performed among all intra-sample p-values
#' to call a mutation.
#' For users wishing to develop their own analysis for other sequencing panel, it requires recalibrated BAM files control samples to be processed to compute the
#' Position-Error Rates stored in a file specified in \code{ber.ctrl.file}.
#'
#' @param patient.dir, char, foldername containing the rBAM folder of the patients. The typical folder hierarchy will consist of 'Plasma/rBAM'
#' @param patient.name, char, filename of the patient .bam file(s) (default NULL read all patients in folder \code{patient.dir})
#' @param pos_ranges.file, char, name of the Rdata file containing the three variables \code{pos_ind}, \code{pos_snp}, \code{pos_ranges} as build by the function \code{PrepareLibrary}. Default NULL, use the position_ranges.rda provides that we used for our analysis.
#' @param ber.ctrl.file, char, pathname of the file providing the background error rates obtained from the controls (default NULL use the provided background error rates obtained from our 29 controls). See \code{background_error_rate.txt} data and \code{BuildCtrlErrorRate} function.
#' @param bai.ext, char, filename extension of the bai files (default '.bai')
#' @param n.trim, integer, number of base positions trimmed at the ends of each amplicon (default 8)
#' @param cov.min, integer, minimal coverture required at each position (default 0)
#' @param force, boolean, (default FALSE) if TRUE force all computations to all files including already processed ones
#' @param show.more, boolean, (default FALSE show only detected positions) if TRUE additional annotations on result plots are given for non-significant mutations
#' @param qcutoff.snv, numeric, proportion of kept base positions ranged by increasing 95th percentile SNV PER in control samples (default 0.95)
#' @param qcutoff.indel, numeric, proportion of kept base positions ranged by increasing 95th percentile INDEL PER in control samples (default 0.99)
#' @param cutoff.sb.ref, numeric, exclude reference positions without cutoff < strand bias < 1-cutoff (default 0.1) (see Supplementary Materials in References)
#' @param cutoff.sb.hotspot, numeric, exclude hotspot positions with Symmetric Odds Ratio test > cutoff (default 3.1) (see Supplementary Materials in References)
#' @param cutoff.sb.nonhotspot, numeric, exclude non-hotspot positions with Symmetric Odds Ratio test > cutoff (default 2.5) (see Supplementary Materials in References)
#' @param hotspot.indel, char, a vector containing the known positions of hotspot deletion/insertion defined as chrX:start:end (default 'chr7:55227950:55249171')
#' @param output.dir, char, name of the folder to save results  (default \code{patient.dir}).
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
#'
#' @return the number of processed patients
#'
#' @export DetectPlasmaMutation
#'
#' @examples
#'      patient.dir=system.file("extdata","4test_only/case/",package="PlasmaMutationDetector")
#'      if (substr(patient.dir,nchar(patient.dir),nchar(patient.dir))!='/')
#'        patient.dir = paste0(patient.dir,'/') # TO RUN UNDER WINDOWS
#'      DetectPlasmaMutation(patient.dir,output.dir=paste0(tempdir(),'/'))
#'
#'
DetectPlasmaMutation = function(patient.dir='./',patient.name=NULL,pos_ranges.file=NULL,ber.ctrl.file=NULL,bai.ext='.bai',
                               n.trim=8,cov.min=0,force=FALSE,show.more=FALSE,
                               qcutoff.snv=0.95,qcutoff.indel=0.99,
                               cutoff.sb.ref=0.1,cutoff.sb.hotspot=3.1,cutoff.sb.nonhotspot=2.5,
                               hotspot.indel='chr7:55227950:55249171',output.dir=patient.dir) {

  ### read and adapt background error rate
  background_error_rate = LoadBackgroundErrorRate(pos_ranges.file,ber.ctrl.file,n.trim)

  ### Identify the so-called noisy SNV and INDEL genomic positions having background error rate higher than an upper quantile.
  ### Quantiles are defined by the cut-off qcutoff.snv and qcutoff.indel
  noisy.snv = background_error_rate[UpperPart(background_error_rate$up.sain, qcutoff.snv),'chrpos']
  noisy.indel = background_error_rate[UpperPart(background_error_rate$indel.up.sain,qcutoff.indel),'chrpos']

  ### extract counts from BAM files
  patient.files = MAF_from_BAM(study.dir=patient.dir,input.filenames=patient.name,bai.ext=bai.ext,
                               pos_ranges.file=pos_ranges.file,force=force,output.dir=output.dir)

  ### directory for the results
  patient.results = paste0(output.dir,'Results/')
  if (!dir.exists(patient.results)) dir.create(patient.results)

  for (i in 1:length(patient.files)) {

    print(paste('Detecting mutations from', patient.files[i]))

    patient.name = rev(strsplit(patient.files[i],'/')[[1]])[1]
    patient.pdf = gsub('_MAF.txt$','.pdf',paste0(patient.results,patient.name))
    patient.hit = gsub('_MAF.txt$','_infos.txt',paste0(patient.results,patient.name))
    patient.mut = gsub('_MAF.txt$','_MUT.txt',paste0(patient.results,patient.name))

    if ((!force & file.exists(patient.hit)) && (file.size(patient.hit)>0)) next() # already processed go to next

    ### Load patient file
    sample = read.table(patient.files[i], header=T, sep="") ### read sample counts
    sample = sample[which(!is.na(sample$MA)),] ### remove positions which haven't been read
    depth = mean(sample$cov,na.rm=T)

    ### Find colums of interest
    ATGC.ind = which(colnames(sample) %in% c('A','T','G','C'))
    ATGC = colnames(sample)[ATGC.ind]
    INDEL.ind = which(colnames(sample) %in% c('longdel','longins'))
    INDEL = colnames(sample)[INDEL.ind]

    ### merge control and sample counts
    samplenoise = merge(sample, background_error_rate, by="chrpos", all=FALSE)

    ### local function made out of samplenoise
    FindColumn = function(i,cols.ind,col.ref,col.names=NULL){
      ind = which(samplenoise[i,cols.ind]==samplenoise[i,col.ref])[1]
      if (is.null(col.names)) {
        ind = cols.ind[ind]
      } else{
        ind = col.names[ind]
      }
      ind
    }

    ### Compute p-values
    samplenoise$probSNV = apply(cbind(c2n(samplenoise$MA),c2n(samplenoise$cov),c2n(samplenoise$E0),c2n(samplenoise$N0)),1, pvalue) # for SNV
    samplenoise$probINDEL = apply(cbind(c2n(samplenoise$longINDEL),c2n(samplenoise$cov),c2n(samplenoise$E0indel),c2n(samplenoise$N0)),1, pvalue) # for INDEL

    ### At each position, start by using SNV reference as default values:
    samplenoise[,'type'] = 'SNV'
    samplenoise[,'hprob'] = samplenoise[,'probSNV']
    ### Identify the REFERENCE allele and the mutated allele:
    for (q in 1:nrow(samplenoise)) {
      samplenoise$Baseref[q] = FindColumn(q,ATGC.ind,'maj',ATGC)
      samplenoise$Basemut[q] = FindColumn(q,ATGC.ind,'MA',ATGC)
      # ### 1.7.0 : to deal with tumoral fraction higher than 50%
      if (as.character(samplenoise$Basemut[q])==as.character(samplenoise$REF[q])) {
        samplenoise$Basemut[q] = as.character(samplenoise$Baseref[q])
        samplenoise$Baseref[q] = as.character(samplenoise$REF[q])
        samplenoise$probSNV[q] = pvalue(c2n(samplenoise$maf[q]),c2n(samplenoise$cov[q]),c2n(samplenoise$E0[q]),c2n(samplenoise$N0[q]))
      }
    }

    ### Correct the INDEL that is the non SNV
    is.INDEL = (samplenoise$probINDEL < samplenoise$probSNV)
    ### Decipher the most significant alteration between SNV, DELETION>2bp or INSERTION>2bp
    is.longins = (samplenoise$longins > samplenoise$longdel)
    samplenoise[which(is.INDEL & is.longins),'type'] = "insertion"
    samplenoise[which(is.INDEL & !is.longins),'type'] = "deletion"
    samplenoise[which(is.INDEL),'hprob'] = samplenoise[which(is.INDEL),'probINDEL']

    ### At each position, compute Strand Bias scores for ref, major variant and indel
    for (q in 1:nrow(samplenoise)) {
      ### compute reference strand bias
      samplenoise[q,'SBref'] = samplenoise[q,paste0(samplenoise[q,'Baseref'],"_SB")]
      refplus = samplenoise[q,paste0(samplenoise[q,'Baseref'],'_plus')]
      refmin = samplenoise[q,paste0(samplenoise[q,'Baseref'],'_minus')]

      ### compute strand bias scores ...
      ### ... for major alteration
      if (samplenoise[q,'type']=='SNV') {
        altplus = samplenoise[q,paste0(samplenoise[q,'Basemut'],'_plus')]
        altmin = samplenoise[q,paste0(samplenoise[q,'Basemut'],'_minus')]
        samplenoise[q,'SBscoreSNV'] = sbscore(refplus, refmin, altplus, altmin)
        samplenoise[q,'hSB'] = samplenoise[q,paste0(samplenoise[q,'Basemut'],"_SB")] # WAS samplenoise[q,'SBsnv']
      }

      ### ... for INDEL
      if (samplenoise[q,'type']!='SNV') {
        altplus = samplenoise[q,paste0(samplenoise[q,'type'],'_plus')]
        altmin = samplenoise[q,paste0(samplenoise[q,'type'],'_minus')]
        samplenoise[q,'SBscoreINDEL'] = sbscore(refplus, refmin, altplus, altmin)
        samplenoise[q,'hSB'] = samplenoise[q,paste0(FindColumn(q,INDEL.ind,'longINDEL',INDEL),"_SB")] # WAS samplenoise[q,'SBindel']
      }

    }

    ### Correction for graphical outputs
    ### Correct the insertions
    samplenoise[which(samplenoise$type=='insertion'),'Baseref'] = 'INS'
    samplenoise[which(samplenoise$type=='insertion'),'Basemut'] = '>2N'

    ### Correct the deletions: length and strand biais
    ### Need to gather consecutive deletions>2bp (sharing the same contig number) into one unique large deletion table
    ind.deletion = which(samplenoise$type=='deletion')
    if (length(ind.deletion)>0) {
      deletion = samplenoise[ind.deletion,]
      for (z in unique(deletion$contig)) {
        deletion$Baseref = 'DEL'
        bool = (deletion$contig==z)
        ind = which(bool)
        if (length(ind)>0) {
          deletion[ind,'Basemut'] = c2n(max(deletion[ind,][[3]])) + 1 - c2n(min(deletion[ind,][[3]])) # gives the number of deleted bp
          deletion[ind,'hprob'] = min(deletion[ind,'probINDEL']) # retains the most significant p-value (p.INDEL) on the large deletion
          deletion[ind,'SBscoreINDEL'] = min(deletion[ind,'SBscoreINDEL']) # retains the lowest Symmetric Odds Ratio Test on the large deletion
          deletion[which(bool & deletion$chrpos %in% noisy.indel),'hprob'] = NA
        }
      }
      deletion = deletion[which(deletion$Basemut!=1),] #ensure that no deletion is <3bp
      samplenoise = rbind(samplenoise[-which(samplenoise$type=='deletion'),], deletion[!duplicated(deletion$contig),])
    }

    ### Remove positions with reference showing outlier strand bias
    samplenoise[which(samplenoise$SBref > 1-cutoff.sb.ref | samplenoise$SBref < cutoff.sb.ref), 'hprob'] = NA

    ### Remove positions showing too high strand bias score (Symmetric Odds Ratio Test)
    samplenoise[which(samplenoise$SBscoreINDEL > cutoff.sb.hotspot), 'hprob'] = NA # for INDEL
    samplenoise[which((samplenoise$hotspot=="Hotspot") & (samplenoise$SBscoreSNV > cutoff.sb.hotspot)), 'hprob'] = NA # for SNV in hotspot
    samplenoise[which((samplenoise$hotspot=="Non-hotspot") & (samplenoise$SBscoreSNV > cutoff.sb.nonhotspot)), 'hprob'] = NA # for SNV in non-hotspot

    ### Remove positions showing 5% highest background error rate in SNV or in INDEL
    type.snv = (samplenoise$type == 'SNV')
    samplenoise[which(type.snv & samplenoise$chrpos %in% noisy.snv ), 'hprob'] = NA
    samplenoise[which(!type.snv & samplenoise$chrpos %in% noisy.indel), 'hprob'] = NA

    ### keep non removed positions (hprob not NA)
    samplenoise = samplenoise[-which(is.na(samplenoise$hprob)),]
    type.snv = (samplenoise$type == 'SNV')

    ### compute the frequency of the mutation
    samplenoise[type.snv,'Freq'] = samplenoise[type.snv,'MA']/samplenoise[type.snv,'cov']
    samplenoise[!type.snv,'Freq'] = samplenoise[!type.snv,'longINDEL']/samplenoise[!type.snv,'cov']

    ### First outlier detection
    samplenoise$outlier="no"
    samplenoise[which(samplenoise$chrpos %in% find.outliers(samplenoise$hprob, samplenoise$chrpos)),'outlier']="outlier"
    nb.outliers=length(samplenoise$outlier[which(samplenoise$outlier!="no")])

    ### Second, we deal now with Hotspots on the smallest p-values
    # find the 250 smallest p-values
    ind = order(samplenoise$hprob, decreasing=FALSE)[1:min(250,nrow(samplenoise))]
    small.pval = samplenoise[ind,'hprob']
    log.small.pval = minus.logit(small.pval)
    #
    # find the upper outliers of a robust boxplot of the finite logit p-values
    upper.wisker = robustbase::adjbox(log.small.pval[is.finite(log.small.pval)], plot=F)$stats[5,1]
    uppers = ind[log.small.pval > upper.wisker]
    for (q in uppers) {
      was.outlier = (samplenoise[q,'outlier']!='no')
      if (samplenoise$type[q]=='SNV' & samplenoise$hotspot[q]=='Hotspot') samplenoise[q,'outlier']="outlier_hotspot"
      if (samplenoise$type[q]!='SNV')
        for (pos in hotspot.indel) { # check if the position is a known indel hotspot positions
          tmp = regmatches(pos,gregexpr(':',pos),invert=T)[[1]]
          if (samplenoise$seqnames[q]==tmp[1] & c2n(samplenoise$start[q])>=c2n(tmp[2]) & c2n(samplenoise$start[q])<=c2n(tmp[3]))
            samplenoise[q,'outlier']="outlier_hotspot"
        }
      if (!was.outlier & samplenoise[q,'outlier']=="outlier_hotspot") nb.outliers=nb.outliers+1
    }

    write.table(samplenoise, file=patient.hit, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    mutations = samplenoise[which(samplenoise$outlier!='no'),c('chrpos','type','Baseref','Basemut','Freq','hotspot','hprob')]
    colnames(mutations) = c('Genomic pos','Mutation type','Base ref','Base mut','Allelic freq','hotspot','p-value')
    write.table(mutations, file=patient.mut, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

    ### For plotting purpose: all possible mutation types are added for legend consistency
    rownames(samplenoise) = 1:nrow(samplenoise)
    samplenoise[nrow(samplenoise)+1,c('hotspot','outlier','type')]=c('Hotspot','outlier','deletion')
    samplenoise[nrow(samplenoise)+1,c('hotspot','outlier','type')]=c('Non-hotspot','no','insertion')
    samplenoise[nrow(samplenoise)+1,c('hotspot','outlier','type')]=c('Non-hotspot','outlier_hotspot','insertion')

    ### For plotting pupose: make type and outlier as.factor with ordered levels
    USNT = unique(samplenoise$type) # YVES 5/09/2016
    samplenoise$type = factor(samplenoise$type, levels = USNT[order(USNT, decreasing=TRUE)]) # YVES 5/09/2016
    USNO = unique(samplenoise$outlier) # YVES 5/09/2016
    samplenoise$outlier = factor(samplenoise$outlier, levels = USNO[order(USNO)]) # YVES 5/09/2016

    ### Plot generation using ggplot2
    grob = grobTree(textGrob(paste0("Mean Depth of Coverage: ", floor(depth),"x"),
                             x=0.45,  y=0.97, hjust=0, gp=gpar(col="red", fontsize=9, fontface="italic")))

    show.pos = (samplenoise$outlier=="outlier_hotspot" | samplenoise$outlier=="outlier")
    if (show.more)
      show.pos = show.pos | ((samplenoise$hotspot=="Hotspot" | samplenoise$Baseref=="DEL" | samplenoise$Baseref=="INS") & samplenoise$hprob<0.000001 & samplenoise$hSB <= 1-cutoff.sb.ref & samplenoise$hSB >= cutoff.sb.ref)

    ind = rev(gregexpr('/',patient.name)[[1]])[1]
    Freq4plot = substr(paste(samplenoise$Freq), 0, 6) # keep 3 decimal for the frequency of the mutation
    pdf = ggplot(samplenoise, aes(sapply(samplenoise$hprob, function(x) -log(max(1e-300,x))), samplenoise$hSB)) +
      geom_point(aes(colour=samplenoise$hotspot, size=samplenoise$outlier, shape=samplenoise$type), na.rm=T) +
      geom_text(aes(label=ifelse(show.pos,
            paste0(as.character(samplenoise$chrpos),".",samplenoise$Baseref,">",samplenoise$Basemut, " AF ", Freq4plot),
            '')),hjust=0, vjust=0.5, angle=65, size=2.5, na.rm=T)  +
      ylab("Strand Biais") + xlab(paste0("-log(p-value)")) + theme_set(theme_bw(base_size = 10)) +
      labs(title = paste(nb.outliers,"Outlier(s)", "\n", substr(patient.name, ind+1, ind+26))) + scale_y_continuous(limits=c(0, 1.0))  +
      scale_x_continuous(limits=c(0, 10e2) ) + annotation_custom(grob) + scale_shape_manual(name="Class", values = c(16, 25, 17)) +
      scale_colour_brewer(name="Hotspot", palette="Set1") + scale_size_manual(name="Outlier", values=c(0.8,5,3))

    # Export results plot and table
    ggsave(pdf, filename=patient.pdf, width = 10, height = 10)

  } # end of loop over patients

  length(patient.files)

}


#' The package provide the SNV and INDEL PERs computed for the Ion AmpliSeq™ Colon and Lung Cancer Panel v2 from 29 controls in a table available in the data file \code{background_error_rate.txt}.
#'
#' This table contains 9 variables for each genomic position
#' \itemize{
#'   \item \code{chrpos}, char, of the form chrN:XXXXXXXXX defining genomic position
#'   \item \code{N0}, integer, the coverture in the controls
#'   \item \code{E0}, integer, the number of errors in the controls
#'   \item \code{p.sain}, numeric,  the ratio E0/N0
#'   \item \code{up.sain}, numeric, the 95th quantile of the Binomial with parameter N0 and E0/N0
#'   \item \code{E0indel}, integer, the amount of indel
#'   \item \code{indel.p.sain}, numeric, the ration E0indel/N0
#'   \item \code{indel.up.sain}, numeric, the 95th quantile of the Binomial with parameter N0 and E0indel/N0
#'   \item \code{hotspot}, char, either 'Non-hotspot' or 'Hotspot' depending if the genomic position is known as hotspot or not.
#' }
#
#' @docType data
#' @usage data(background_error_rate)
#' @name background_error_rate
#' @aliases background_error_rate.txt
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @keywords data
#' @seealso \code{BuildCtrlErrorRate}
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
NULL

#' The package provide the positions and ranges computed for the Ion AmpliSeq™ Colon and Lung Cancer Panel v2 as a Rdata file \code{positions_ranges.rda}.
#'
#' This file contains 4 variables
#' \itemize{
#'   \item \code{pos_ind}, vector of chars, of the form chrN:XXXXXXXXX defining genomic positions of the Ion AmpliSeq™ Colon and Lung Cancer Panel v2
#'   \item \code{pos_snp}, vector of chars, of the form chrN:XXXXXXXXX defining the known snp genomic positions
#'   \item \code{pos_ranges}, GRanges object, describing the 92 amplicons of the Ion AmpliSeq™ Colon and Lung Cancer Panel v2
#' }
#
#' @docType data
#' @usage data(positions_ranges)
#' @name positions_ranges
#' @aliases positions_ranges.rda pos_ind pos_snp pos_ranges
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @keywords data
#' @seealso \code{Prepare_Library}
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
NULL

#' The package provide a list of known hotspot positions located on the amplicons of the Ion AmpliSeq™ Colon and Lung Cancer Panel v2 as a txt file \code{hotspot.txt} which contains a vector/variable ---named chrpos (first row)--- of chars, of the form chrN:XXXXXXXXX defining genomic positions.
#'
#
#' @docType data
#' @usage data(hotspot)
#' @name hotspot
#' @aliases hotspot.txt
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @keywords data
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
NULL

