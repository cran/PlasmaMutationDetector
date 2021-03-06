% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{MAF_from_BAM}
\alias{MAF_from_BAM}
\title{function MAF_from_BAM}
\usage{
MAF_from_BAM(study.dir = "Plasma/", input.filenames = NULL,
  bai.ext = ".bai", pos_ranges.file = NULL, force = FALSE,
  output.dir = study.dir)
}
\arguments{
\item{study.dir, }{char, name of the folder containing the rBAM directory  (default 'Plasma/'). The typical folder hierarchy will consist of 'Plasma/rBAM'}

\item{input.filenames, }{a vector of char (default NULL), the names of the BAM files to process. If NULL all BAM files in the rBAM folder will be processed}

\item{bai.ext, }{char, filename extension of the bai files (default '.bai')}

\item{pos_ranges.file, }{char, name of the Rdata file containing the three variables \code{pos_ind}, \code{pos_snp} and \code{pos_ranges} as build by the function \code{PrepareLibrary}. Default NULL, use the position_ranges.rda provided, used for our analysis.}

\item{force, }{boolean, (default FALSE) if TRUE force all computations to all files including already processed ones}

\item{output.dir, }{char, name of the folder to save results  (default \code{study.dir}).}
}
\value{
the path/names of the MAF files
}
\description{
Read BAM files  and create MAF file. BAMfiles are stored in a sub-folder '/rBAM'.
MAF files are intermediate files stored in a sub-folder '/BER'.
MAF files contain the raw counts of A,T,C,G, insertion, deletion, insertion>2bp, deletion >2bp for strand plus and stand minus.
Note : we strongly recommand to externally recalibrate BAM files using tools like GATK.
}
\examples{
  \dontrun{
     ctrl.dir = system.file("extdata", "4test_only/ctrl/", package = "PlasmaMutationDetector")
     if (substr(ctrl.dir,nchar(ctrl.dir),nchar(ctrl.dir))!='/')
       ctrl.dir = paste0(ctrl.dir,'/') # TO RUN UNDER WINDOWS
     MAF_from_BAM(ctrl.dir,force=TRUE,output.dir=paste0(tempdir(),'/'))
   }

}
\references{
\emph{Analysis of base-position error rate of next-generation sequencing to detect tumor mutations in circulating DNA} N. Pécuchet, Y. Rozenholc, E. Zonta, D. Pietraz, A. Didelot, P. Combe, L. Gibault, J-B. Bachet, V. Taly, E. Fabre, H. Blons, P. Laurent-Puig in \emph{Clinical Chemistry}
}
\author{
N. Pécuchet, P. Laurent-Puig and Y. Rozenholc
}
