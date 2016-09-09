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
#' @importFrom stats dnorm integrate median pbinom quantile
#' @importFrom utils data download.file read.delim read.table write.table
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

minus.logit = function(P) sapply(P,function(x) {if(is.nan(x)||(x<1)) return(-log(x/(1-x))) else return(Inf)})

pvalue = function(E,N=NULL,E0=NULL,N0=NULL,with.info=FALSE) {
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
  integrate(psi,lower=-30,upper=30,subdivisions=10000L)$value
}


find.outliers = function(P,chr.pos,thr.R=6.5,thr.P=0.000001,with.info=FALSE) {

  logit = function(P) sapply(P,function(x) if(is.nan(x)||(x<1)) return(log(x/(1-x))) else return(Inf))

  P[which(P==1)]=0.999999999
  nP = length(P)
  names(P) = chr.pos
  L.tri = sort(-logit(P))
  L.tri[L.tri==Inf] = 10*max(L.tri[L.tri!=Inf])
  dL = diff(L.tri);
  L.tri.short=L.tri[-1];
  L.moy=L.tri[-nP];
  L.moy=sign(L.moy)*(abs(L.moy))^0.5
  Ratio = dL/L.moy
  outliers = which((Ratio>thr.R)&(L.tri.short>-logit(thr.P)))

  ### ADD YVES 29/03/2016
  if (length(outliers)>0) outliers = which(L.tri.short>=min(L.tri.short[min(outliers)]))

  ind = names(outliers)
  if (with.info) {
    plot(L.tri[-1][L.moy>0],Ratio[L.moy>0],col=1)
    abline(v=-logit(thr.P),col='red'); abline(h=thr.R,col='green'); abline(h=6,col='cyan'); abline(h=7,col='violet');
    if (length(outliers)>0) outliers = which(L.tri.short>=min(L.tri.short[min(outliers)]))
    if (length(outliers)>0) points(L.tri[-1][outliers],Ratio[outliers],col='red',pch='x')
    if (length(outliers)>0) text(L.tri[-1][outliers],Ratio[outliers],labels=ind,cex=0.6,pos=2,col=1)
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
#' @param output.name, char, filename to save \code{pos_ind} and \code{pos_snp} in \code{info.dir} (default 'positions_ranges.rda')
#' @param load.from.broad.insitute, boolean, if TRUE load \code{snp.filename} from Broad Institute ftp server otherwise use the file positions_ranges_broad.rda (default FALSE)
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @seealso positions_ranges,
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
#'
#' @return Save the following variables in a .rda file defined by \code{output.name} in the folder defined by \code{info.dir}:
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
#'    PrepareLibrary(info.dir='./',snp.extra=bad.pos)
#'
#'
PrepareLibrary = function(info.dir='Info/',
                          bed.filename='lungcolonV2.bed.txt', # amplicon positions
                          snp.filename='ExAC.r0.3.sites.vep.vcf.gz', # known snp
                          snp.extra=c("chr4:1807909","chr7:140481511", "chr18:48586344", "chr14:105246474", "chr19:1223055"), # extra snp
                          output.name='positions_ranges.rda',
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


  output.name = paste0(info.dir,output.name)
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
#' @param pos_ranges.file, char, name of the Rdata file containing the three variables \code{pos_ind}, \code{pos_snp} and \code{pos_ranges} as build by the function \code{PrepareLibrary}. Default NULL, use the position_ranges.rda provided, used for our analysis.
#' @param bai.ext, char, filename extension of the bai files (default '.bai')
#' @param force, boolean, (default FALSE) if TRUE force all computations to all files including already processed ones
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
#'
#' @return the number of processed bam files
#'
#' @export MAF_from_BAM
#'
#' @examples
#'    \dontrun{
#'      ctrl.dir = system.file("extdata", "4test_only/ctrl/", package = "PlasmaMutationDetector")
#'      if (substr(ctrl.dir,nchar(ctrl.dir),nchar(ctrl.dir))!='/')
#'        ctrl.dir = paste0(ctrl.dir,'/') # TO RUN UNDER WINDOWS
#'      MAF_from_BAM(ctrl.dir,force=TRUE)
#'    }
#'
MAF_from_BAM = function(study.dir='Plasma/',input.filenames=NULL,bai.ext='.bai',
                        pos_ranges.file=NULL,force=FALSE) {

  pos_ind = pos_snp = pos_ranges = NULL

  if (is.null(pos_ranges.file)) {
    data(positions_ranges,envir=environment())
  } else {
    if (!file.exists(pos_ranges.file)) stop(paste('STOP :::',pos_ranges.file,' does not exist'))
    load(pos_ranges.file)
  }

  ######## Input/Output folders
  bam.dir = paste0(study.dir,'rBAM/')        # Folder with BAM files and BAI index files
  ber.dir = paste0(study.dir,'BER/')         # Folder for Background Error Rate files extracted from BAM files
  if (!dir.exists(ber.dir)) dir.create(ber.dir)

  ######## pileup parameters
  pileup.param = PileupParam(max_depth=100000, min_base_quality=20, min_mapq=5, min_nucleotide_depth=1, min_minor_allele_depth=0,
                             distinguish_strands=TRUE, distinguish_nucleotides=TRUE, ignore_query_Ns=TRUE, include_deletions=TRUE,
                             include_insertions=TRUE)

  ######## Samples to process
  if (is.null(input.filenames)) {
    print(paste0("Looking for bam and bai in ",bam.dir,' ...'))
    input.filenames = dir(bam.dir, pattern="*.bam$") # list all BAM files in bam.dir
  }

  bam.files = paste0(bam.dir,input.filenames) # build bam.files as paste0(study.dir,'rBAM/',input.filenames[i])

  if (length(bam.files)==0) stop(paste('STOP ::: no bam file found ... NOTE: bam files need to be in a subdirectory called rBAM'))

  for( g in 1:length(bam.files)){

    bam.name = bam.files[g]
    bai.name = gsub('.bam$',bai.ext,bam.name)
    ber.name = gsub(bam.dir,ber.dir,gsub('.bam$','_MAF.txt',bam.name),fixed=TRUE)
    berzip.name = paste0(ber.name,'.zip')

    if (!force & file.exists(ber.name)) next()
    if (!force & file.exists(berzip.name)) {
      write.table(readLines(unz(berzip.name,rev(strsplit(ber.name,'/')[[1]])[1])),file=ber.name,sep=' ',row.names=FALSE,col.names=FALSE,quote=FALSE)
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
    ft = scanBam(bf, scanBamParam=ScanBamParam(which=pos_ranges, what=scanBamWhat()),pileupParam=pileup.param)

    ### Count SNV from the BAM file
    res = pileup(bf, scanBamParam=ScanBamParam(which=pos_ranges, what=scanBamWhat()),pileupParam=pileup.param)

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

      ### assign the same contig id to deleted positions if they are distant of 3bp and consistently frequent
      for (i in 1:(nrow(contig)-2))
        if ( (as.numeric(contig[i,'start']) + 2) == as.numeric(contig[i+2,'start']) )
          if(abs(contig[i,'-'] - contig[i+2,'-'])/max(contig[i,'-'], contig[i+2,'-']) < 0.2 & abs(contig[i,'-'] - contig[i+1,'-'])/max(contig[i,'-'], contig[i+1,'-']) < 0.2)
            contig[i:i+2,'contig'] = i
      ### remove positions with inconsistent deletions (<2bp or different frequency)
      contig2 = contig[-which(is.na(contig$contig)),]

      if (nrow(contig2)>1)
        for (i in 1:(nrow(contig2)-1))
          for (j in 1:(nrow(contig2)-i))
            if (contig2[i,'contig'] == (contig2[i+j,'contig'] - j) )
              contig2[i+j,'contig']=contig2[i,'contig']

      for (i in 1:(nrow(contig)-1))
        if ( i %in% contig2$contig) {
          contig3 = rbind(contig2[which(contig2$contig == i),][1,],contig2[which(contig2$contig == i),][1,])
          contig3[1,'start']= as.numeric(contig3[1,'start'])-2
          contig3[2,'start']= as.numeric(contig3[2,'start'])-1
          contig2 = rbind(contig2,contig3)
        }

      contig2[,'chrpos']=paste0(contig2[,'seqnames'],":",contig2[,'start'])

      ### A warning to prevent the false duplication of deletions
      if (length(contig2$chrpos[duplicated(contig2$chrpos)])>0)
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
    write.table(select[!duplicated(select$chrpos),], file=ber.name, row.names=FALSE, col.names=TRUE) ####  save temporary files
  }

  return(length(bam.files))

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
#' @author N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
#' @references \emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
#'
#' @return the number of processed files
#'
#' @export BuildCtrlErrorRate
#'
#' @examples
#'    ctrl.dir = system.file("extdata", "4test_only/ctrl/", package = "PlasmaMutationDetector")
#'    if (substr(ctrl.dir,nchar(ctrl.dir),nchar(ctrl.dir))!='/')
#'      ctrl.dir = paste0(ctrl.dir,'/') # TO RUN UNDER WINDOWS
#'    BuildCtrlErrorRate(ctrl.dir)
#'
BuildCtrlErrorRate = function(ctrl.dir='Plasma ctrl/',bai.ext='.bai',pos_ranges.file=NULL,hotspot.file=NULL,
                              force=FALSE) {

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
    hotspot = read.delim(hotspot.file, header=FALSE)
  }

  ### extract counts from all BAM files found in ctrl.dir
  MAF_from_BAM(study.dir=ctrl.dir,input.filenames=NULL,bai.ext=bai.ext,pos_ranges.file=pos_ranges.file,force=force)

  ### list Background-Error-Rate files
  ctrl.files = dir(paste0(ctrl.dir,'BER'), pattern="*.txt$", full.names=TRUE) # list all Background Error Rate files

  ### Assemble background errors from control samples
  X = read.table(ctrl.files[1], header=T, sep=" ")
  if (length(ctrl.files)>1) ind.g = 2:length(ctrl.files) else ind.g = c(1,1) # just for testing purpose with one bam
  for (g in ind.g)
    X = merge(read.table(ctrl.files[g],header=T,sep=" "), X, by="chrpos", all=TRUE)

  ### Few controls to be sure ... useless???
  if (anyDuplicated(X$chrpos)) {
    cat( "!!!!! Warning !!!! Duplicated genomic positions in controls ... remove duplicated positions" )
    X = X[!duplicated(X$chrpos),]
  }

  X = X[(X[,1] %in% pos_ind),] #### ensure to select positions in defined genomic ranges
  X = X[!(X[,1] %in% pos_snp),] #### ensure to remove SNP positions and background error rate for SNPs (E0) and INDELS (E0indel) in samples

  X.col = colnames(X)
  col.cov = grep('cov',X.col) # column indexes of the coverture
  col.del = grep('longdel($|[^_])',X.col) # column indexes of the longdel
  col.ins = grep('longins($|[^_])',X.col) # column indexes of the longins
  col.maj = grep('^maj',X.col) # column indexes of the major

  sum.no.na = function(x) sum(x,na.rm=TRUE)
  q95 = function(x) quantile(x,0.95,na.rm=TRUE)

  N0 = apply(X[,col.cov],1,sum.no.na) # Sum of the depths of all control samples for each position
  E0 = apply(X[,col.cov]-X[,col.maj],1,sum.no.na) # Sum of the SNV background counts of all control samples for each position

  X = X[which(!is.na(N0) | N0 > 5000),]  # Remove genomic positions which are not covered in all control samples or with sum of depths < 5,000x

  ### SNV rates
  p.sain = E0/N0 # Control SNV Position-Error Rate
  up.sain = apply((X[,col.cov]-X[,col.maj])/X[,col.cov],1,q95) # 95th percentile for Control SNV Position-Error Rate

  ### INDEL rates
  E0indel= apply(cbind(X[,col.del], X[,col.ins]),1,sum.no.na) # Sum of the INDEL background counts of all control samples for each position
  indel.p.sain = E0indel/N0 # Control INDEL Position-Error Rate
  indel.up.sain = apply((X[,col.del] + X[,col.ins]) / X[,col.cov], 1, q95) # 95th percentile for Control INDEL Position-Error Rate

  # Joint dataframe containing Control SNV and INDEL Position-Error Rate together with their 95th percentile, and the raw counts data
  chrpos = as.character(X[,'chrpos'])
  background_error_rate <- data.frame(cbind(chrpos,N0,E0,p.sain,up.sain,E0indel,indel.p.sain,indel.up.sain))

  # Precise what are the hotspot positions
  background_error_rate$hotspot <- 'Non-hotspot'
  background_error_rate[which(background_error_rate$chrpos %in% hotspot$V1 ),'hotspot'] = 'Hotspot'

  # Save Background error file...
  write.table(background_error_rate, file=paste0(ctrl.dir,"background_error_rate.txt"), row.names=FALSE, col.names=TRUE)

  return(length(ctrl.files))
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
#' @param force, boolean, (default FALSE) if TRUE force all computations to all files including already processed ones
#' @param show.more, boolean, (default FALSE show only detected positions) if TRUE additional annotations on result plots are given for non-significant mutations
#' @param qcutoff.snv, numeric, proportion of kept base positions ranged by increasing 95th percentile SNV PER in control samples (default 0.95)
#' @param qcutoff.indel, numeric, proportion of kept base positions ranged by increasing 95th percentile INDEL PER in control samples (default 0.99)
#' @param cutoff.sb.hotspot, numeric, exclude variants on hotspot positions with Symmetric Odds Ratio test > cutoff (default 3.1) (see Supplementary Materials in References)
#' @param cutoff.sb.nonhotspot, numeric, exclude variants on non-hotspot positions with Symmetric Odds Ratio test > cutoff (default 2.5) (see Supplementary Materials in References)
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
#'      DetectPlasmaMutation(patient.dir)
#'
#'
DetectPlasmaMutation = function(patient.dir='./',patient.name=NULL,pos_ranges.file=NULL,ber.ctrl.file=NULL,bai.ext='.bai',
                               n.trim=8,force=FALSE,show.more=FALSE,
                               qcutoff.snv=0.95,qcutoff.indel=0.99,cutoff.sb.hotspot=3.1,cutoff.sb.nonhotspot=2.5) {

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

  ####### Identify positions with 95th percentile SNV and INDEL Position-Error Rates (SNV-PER and INDEL-PER)
  ####### larger than cut-off
  noisy.snv = background_error_rate[which(background_error_rate$up.sain > quantile(background_error_rate$up.sain, qcutoff.snv)),'chrpos']
  noisy.indel =
    background_error_rate[which(background_error_rate$indel.up.sain > quantile(background_error_rate$indel.up.sain,qcutoff.indel) ),'chrpos']

  ### extract counts from BAM files
  MAF_from_BAM(study.dir=patient.dir,input.filenames=patient.name,bai.ext=bai.ext,pos_ranges.file=pos_ranges.file,force=force)

  if (is.null(patient.name)) patient.files = dir(paste0(patient.dir,'BER'), pattern="*_MAF.txt$", full.names=TRUE) # list all Background Error Rate files
  else patient.files = paste0(patient.dir,'BER/',gsub('.bam$','_MAF.txt',patient.name))

  patient.results = paste0(patient.dir,'Results/')
  if (!dir.exists(patient.results)) dir.create(patient.results)

  for (i in 1:length(patient.files)) {

    patient.name = patient.files[i]

    print(paste('Detecting mutations from',patient.name))

    patient.pdf = gsub('_MAF.txt$','.pdf',gsub('BER/','Results/',patient.name))
    patient.hit = gsub('_MAF.txt$','_infos.txt',gsub('BER/','Results/',patient.name))
    patient.mut = gsub('_MAF.txt$','_MUT.txt',gsub('BER/','Results/',patient.name))

    if ((!force & file.exists(patient.hit)) && (file.size(patient.hit)>0)) next()

    sample = read.table(patient.name, header=T, sep="") ### read sample counts
    samplenoise = merge(sample, background_error_rate, by="chrpos", all=FALSE) ### merge control and sample counts

    depth = mean(samplenoise$cov,na.rm=T)

    samplenoise = samplenoise[which(!is.na(samplenoise$MA)),] ### remove positions which haven't been read (by mistake)

    ### Compute p-values
    samplenoise$probSNV = apply(cbind(c2n(samplenoise$MA),c2n(samplenoise$cov),c2n(samplenoise$E0),c2n(samplenoise$N0)),1, pvalue) # for SNV
    samplenoise$probINDEL = apply(cbind(c2n(samplenoise$longINDEL),c2n(samplenoise$cov),c2n(samplenoise$E0indel),c2n(samplenoise$N0)),1, pvalue) # for INDEL

    ## At each position, decipher the most significant alteration between SNV, DELETION>2bp or INSERTION>2bp
    ### Use SNV reference as default values:
    samplenoise[,'type'] = 'SNV'
    samplenoise[,'hprob'] = samplenoise[,'probSNV']
    ### Correct the non SNV that is the INDEL:
    is.INDEL = (samplenoise$probINDEL < samplenoise$probSNV)
    is.longins = (samplenoise$longins > samplenoise$longdel)
    samplenoise[which(is.INDEL & is.longins),'type'] = "insertion"
    samplenoise[which(is.INDEL & !is.longins),'type'] = "deletion"
    samplenoise[which(is.INDEL),'hprob'] = samplenoise[which(is.INDEL),'probINDEL']

    ATGC = colnames(samplenoise)[4:7]
    INDEL = colnames(samplenoise)[c(41,46)]

    ## At each position, compute Strand Bias scores
    for (q in 1:nrow(samplenoise)) {
      samplenoise[q,'SBsnv'] = max(samplenoise[q,paste0(sub("\\..*", "\\1",ATGC[samplenoise[q,4:7]==samplenoise[q,'MA']]),"_SB")])
      samplenoise[q,'SBindel'] = max(samplenoise[q,paste0(sub("\\..*", "\\1",INDEL[samplenoise[q,c(41,46)]==samplenoise[q,'longINDEL']]),"_SB")])
      samplenoise[q,'SBref'] = max(samplenoise[q,paste0(sub("\\..*", "\\1",ATGC[samplenoise[q,4:7]==samplenoise[q,'maj']]),"_SB")])

      ### strand bias scores
      ref = ATGC[which.max(samplenoise[q,4:7])]
      refplus = samplenoise[q,paste0(ref,'_plus')]
      refmin = samplenoise[q,paste0(ref,'_minus')]

      # strand bias score for major alteration
      if (samplenoise[q,'type']=='SNV') {
        alt = ATGC[samplenoise[q,4:7]==samplenoise[q,'MA']][1]
        altplus = samplenoise[q,paste0(alt,'_plus')]
        altmin = samplenoise[q,paste0(alt,'_minus')]
        samplenoise[q,'SBscoreSNV'] = sbscore(refplus, refmin, altplus, altmin)
        samplenoise[q,'hSB'] = samplenoise[q,'SBsnv']
      }

      # strand bias score for INDEL
      if (samplenoise[q,'type']!='SNV') {
        altplus = samplenoise[q,paste0(samplenoise[q,'type'],'_plus')]
        altmin = samplenoise[q,paste0(samplenoise[q,'type'],'_minus')]
        samplenoise[q,'SBscoreINDEL'] = sbscore(refplus, refmin, altplus, altmin)
        samplenoise[q,'hSB'] = samplenoise[q,'SBindel']
      }
    }

    #### Used for results plots, identify the REFERENCE allele, the mutated allele if SNV, the deletion length, the insertion, and the strand biais of the alteration
    for (p in 1:nrow(samplenoise)) {
      samplenoise[p,'Baseref']=paste0(sub("\\..*", "\\1",ATGC[samplenoise[p,4:7]==samplenoise[p,'maj']]), collapse='')
      samplenoise[p,'Basemut']=paste0(sub("\\..*", "\\1",ATGC[samplenoise[p,4:7]==samplenoise[p,'MA']]), collapse='')
    }
    samplenoise[which(samplenoise$type=='insertion'),'Baseref'] = 'INS'
    samplenoise[which(samplenoise$type=='insertion'),'Basemut'] = '>2N'

    # Gather consecutive deletions>2bp (sharing the same contig number) into one unique large deletion
    ind.deletion = which(samplenoise$type=='deletion')
    if (length(ind.deletion)>0) {
      deletion = samplenoise[ind.deletion,]
      for (z in unique(deletion$contig)) {
        deletion$Baseref <- 'DEL'
        bool = (deletion$contig==z)
        ind = which(bool)
        if (length(ind)>0) {
          deletion[ind,'Basemut'] <- c2n(max(deletion[ind,][[3]])) + 1 - c2n(min(deletion[ind,][[3]])) # gives the number of deleted bp
          deletion[ind,'hprob'] <- min(deletion[ind,'probINDEL']) # retains the most significant p-value (p.INDEL) on the large deletion
          deletion[ind,'SBscoreINDEL'] <- min(deletion[ind,'SBscoreINDEL']) # retains the lowest Symmetric Odds Ratio Test on the large deletion
          deletion[which(bool & deletion$chrpos %in% noisy.indel),'hprob'] <- NA
        }
      }
      deletion = deletion[which(deletion$Basemut!=1),] #ensure that no deletion is <3bp
      samplenoise = rbind(samplenoise[-which(samplenoise$type=='deletion'),], deletion[!duplicated(deletion$contig),])
    }


    samplenoise[which(samplenoise$SBref > 0.90 | samplenoise$SBref < 0.10), 'hprob'] = NA # Remove positions with extreme strand biais on reference allele
    samplenoise[which(samplenoise$SBscoreSNV > cutoff.sb.hotspot), 'hprob'] = NA # SNV Symmetric Odds Ratio Test cut-off = 3.1
    samplenoise[which(samplenoise$SBscoreINDEL > cutoff.sb.hotspot), 'hprob'] = NA # INDEL Symmetric Odds Ratio Test cut-off = 3.1
    samplenoise[which((samplenoise$hotspot=="Non-hotspot") & (samplenoise$SBscoreSNV > cutoff.sb.nonhotspot)), 'hprob'] = NA # More stringent Symmetric Odds Ratio Test cut-off to remove non-hotspot SNV

    # Remove positions with high 95th percentile Control SNV Position-Error Rate and Control INDEL Position-Error Rate
    type.snv = (samplenoise$type == 'SNV')
    samplenoise[which(type.snv & samplenoise$chrpos %in% noisy.snv ), 'hprob'] = NA
    samplenoise[which(!type.snv & samplenoise$chrpos %in% noisy.indel), 'hprob'] = NA

    samplenoise = samplenoise[-which(is.na(samplenoise$hprob)),]

    type.snv = (samplenoise$type == 'SNV')
    ind = which(type.snv)
    samplenoise[ind,'Freq'] = samplenoise[ind,'MA']/samplenoise[ind,'cov'] # compute the frequency of the mutation
    ind = which(!type.snv)
    samplenoise[ind,'Freq'] = samplenoise[ind,'longINDEL']/samplenoise[ind,'cov'] # compute the frequency of the mutation
    Freq4plot = substr(paste(samplenoise$Freq), 0, 6) # keep 3 decimal for the frequency of the mutation

    outlier = 0
    samplenoise$outlier="no"

    # 1.Hotspots on the smallest p-values
    if (nrow(samplenoise)>10) {
      pboxplot = samplenoise[which(samplenoise$hprob %in% sort(samplenoise$hprob, decreasing=FALSE)[1:250]),]
      pboxplot$logp = sapply(pboxplot$hprob, minus.logit)

      write.table(pboxplot$logp,file = 'pboxplot_hprob_tmp.txt',col.names = F,row.names = F)

      xx = rownames(pboxplot)[which(pboxplot$logp %in% robustbase::adjbox(pboxplot$logp[is.finite(pboxplot$logp)], plot=F)$out == T &
                                      pboxplot$logp > median(pboxplot$logp))]
      xxx = rownames(pboxplot[which(pboxplot$hprob<0.05),])
      xxxx = rownames(pboxplot[which((pboxplot$hotspot=="Hotspot" & pboxplot$type=="SNV") |
                                       (pboxplot$seqnames=='chr7' & c2n(pboxplot$start)>=55227950  &
                                          c2n(pboxplot$start)<=55249171 &
                                          (pboxplot$Baseref=="DEL" | pboxplot$Baseref=="INS"))
      ),]
      )
      samplenoise[xx[which(xx %in% xxx & xx %in% xxxx)],'outlier']="outlier_hotspot"
      outlier=length(xx[which(xx %in% xxx & xx %in% xxxx)])
    }

    # 2.First outlier detection
    samplenoise[which(samplenoise$chrpos %in% find.outliers(samplenoise$hprob, samplenoise$chrpos)),'outlier']="outlier"
    nb.outliers=length(samplenoise$outlier[which(samplenoise$outlier!="no")])

    write.table(samplenoise, file=patient.hit, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    mutations = samplenoise[which(samplenoise$outlier!='no'),c('chrpos','type','Baseref','Basemut','Freq','hotspot','hprob')]
    colnames(mutations) = c('Genomic pos','Mutation type','Base ref','Base mut','Allelic freq','hotspot','p-value')
    write.table(mutations, file=patient.mut, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

    # For plotting purpose, we add all possible types of mutations so that the legend is consistent on all samples, mutated or not
    rownames(samplenoise) = 1:nrow(samplenoise)
#     samplenoise$hotspot = factor(samplenoise$hotspot,levels=c('hotspot','non-hotspot'))
#     samplenoise$outlier = factor(samplenoise$outlier,levels=c('outlier','no','outlier_hotspot'))
#     samplenoise$outlier = factor(samplenoise$type,levels=c('deletion','insertion','snv'))
    samplenoise[nrow(samplenoise)+1,c('hotspot','outlier','type')]=c('Hotspot','outlier','deletion')
    samplenoise[nrow(samplenoise)+1,c('hotspot','outlier','type')]=c('Non-hotspot','no','insertion')
    samplenoise[nrow(samplenoise)+1,c('hotspot','outlier','type')]=c('Non-hotspot','outlier_hotspot','insertion')
    USNT = unique(samplenoise$type) # YVES 5/09/2016
    samplenoise$type = factor(samplenoise$type, levels = USNT[order(USNT, decreasing=TRUE)]) # YVES 5/09/2016
    # samplenoise$type = factor(samplenoise$type, levels = samplenoise$type[order(samplenoise$type, decreasing=TRUE)])
    USNO = unique(samplenoise$outlier) # YVES 5/09/2016
    samplenoise$outlier = factor(samplenoise$outlier, levels = USNO[order(USNO)]) # YVES 5/09/2016
    # samplenoise$outlier = factor(samplenoise$outlier, levels = samplenoise$outlier[order(samplenoise$outlier)])

    # The final plot is generated using package ggplot2
    grob = grobTree(textGrob(paste0("Mean Depth of Coverage: ", floor(depth),"x"),
                             x=0.45,  y=0.97, hjust=0, gp=gpar(col="red", fontsize=9, fontface="italic")))

    show.pos = (samplenoise$outlier=="outlier_hotspot" | samplenoise$outlier=="outlier")
    if (show.more)
      show.pos = show.pos | ((samplenoise$hotspot=="Hotspot" | samplenoise$Baseref=="DEL" | samplenoise$Baseref=="INS") & samplenoise$hprob<0.00001 & samplenoise$hSB < 0.75 & samplenoise$hSB > 0.25)

    ind = rev(gregexpr('/',patient.name)[[1]])[1]
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

    samplenoise = NULL

  }

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

#' The package provide a list of known hotspot positions located on the amplicons of the Ion AmpliSeq™ Colon and Lung Cancer Panel v2 as a txt file \code{hotspot.txt} which contains a vector of vector of chars, of the form chrN:XXXXXXXXX defining genomic positions.
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

