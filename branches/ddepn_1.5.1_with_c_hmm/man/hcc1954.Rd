\name{hcc1954}
\Rdversion{1.1}
\alias{hcc1954}
\alias{hcc1954raw}
\alias{hcc1954fcf}
\alias{hcc1954fcfmedian}
\docType{data}
\title{
hcc1954
}
\description{
Example phosphoproteomic data set. Phosphorylation of 16 proteins from the ERBB signalling
network were measured in 10 time points (between 0 and 60 minutes) in HCC1954 cell line.
Measurements were generated via Reverse Phase Protein Array technique.
}
\format{Matrix 16x540. Rows hold proteins, Columns hold the Experiments and time points.\cr
\emph{Column naming:} STIMULUS_biolreplicate_time\cr
\emph{Rows naming:} Protein_Phosphorylation site.\cr
}
\details{
 Use \code{data(hcc1954)} to obtain the FCF normalised and replicate centered data matrix or
 \code{data(hcc1954raw)} to obtain the raw data matrices. The latter comprise the 
 raw data \code{hcc1954raw}, the FCF values \code{hcc1954fcf} and the medians of the 
 experiment FCF values \code{hcc1954fcfmedian}. See examples for how to transform the raw 
 to the normalised data.\cr
 
  
 The following experiments were performed: stimulation with ...\cr 
  ...EGF\cr
  ...HRG\cr
  ...EGF&HRG\cr

 The following proteins were measured:
 
 \tabular{lc}{
    \emph{protein_phosphorylationsite} \tab entrez id \cr
    pEGFR_Y1068 \tab 1956\cr
    pERBB2_Y1112  \tab  2064\cr
    pERK12_T202Y204 \tab 5595  5594\cr
    pAKT_S473   \tab   207  208   \cr
    pPDK1_S241 \tab 5170\cr
    pMEK_S217S221  \tab  5604 5605\cr 
    pPLCgamma_S1248 \tab 5335  5336\cr
    pPKCalpha_S657Y658 \tab 5578\cr
    pp38_T180Y182 \tab 1432\cr
    pSRC_Y416   \tab  6714\cr
    pmTOR_S2448 \tab 2475\cr
    pp70S6K_T389   \tab  6198\cr 
    pGSK3_Y279Y216 \tab 2931 2932\cr
    pPRAS_T246   \tab     100136753\cr
    pERBB3_Y1289 \tab2065\cr
    pERBB4_Y1162 \tab 2066
 }
}

\note{
To use these data matrices as inputs for \code{\link{ddepn}}, the information on the 
biological replicates in the column names have to be removed using the function
\code{\link{format_ddepn}}, since all replicates are treated equally during inference. 
The information on the biological replicates was left in the data matrix, however, 
since it might be useful for different approaches, too.
}

\seealso{
\code{\link{format_ddepn}}
\code{\link{center_ddepn}}
}

\references{
Bender et. al. 2010: Dynamic deterministic effects propagation networks: learning 
 signalling pathways from longitudinal protein array data; Bioinformatics,
 Vol. 26(18), pp. i596-i602
}
\examples{
\dontrun{

library(ddepn)

## get the normalised data matrix
data(hcc1954)
colnames(hcc1954)
rownames(hcc1954)

## get the raw matrices
data(hcc1954raw)

## perform FCF normalisation
datfcf <- hcc1954raw / hcc1954fcf * hcc1954fcfmedian

## center the data, datc is then equal to the hcc1954 matrix
datc <- center_ddepn(datfcf)
}
}
\keyword{datasets}
