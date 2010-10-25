\name{kegggraphs}
\Rdversion{1.1}
\alias{kegggraphs}
\docType{data}
\title{
kegggraphs
}
\description{
Contains a list of selected KEGG pathways, downloaded in October 2010.
}
\format{List of 78, each Element holding the KEGG pathway id, a graphNEL
 object and an adjacency matrix.}
\details{
The KEGG Hierarchy was downloaded and transformed into graphNEL objects
by using KEGGgraph. The helper \code{\link{kegggraph.to.detailed.adjacency}} 
can be used to convert the graphNEL objects into simple adjacency matrices.
A selection of pathways was performed to restrict to the signalling and disease
pathways:\cr

The following pathways are included:\cr
\tabular{ll}{
KEGG pathway ID \tab Name\cr
04010 \tab MAPK signaling pathway \cr
04012 \tab ErbB signaling pathway \cr
04020 \tab Calcium signaling pathway \cr
04060 \tab Cytokine-cytokine receptor interaction \cr
04080 \tab Neuroactive ligand-receptor interaction \cr
04110 \tab Cell cycle \cr
04115 \tab p53 signaling pathway \cr
04120 \tab Ubiquitin mediated proteolysis \cr
04130 \tab SNARE interactions in vesicular transport \cr
04140 \tab Regulation of autophagy \cr
04150 \tab mTOR signaling pathway \cr
04210 \tab Apoptosis \cr
04310 \tab Wnt signaling pathway \cr
04320 \tab Dorso-ventral axis formation \cr
04330 \tab Notch signaling pathway \cr
04340 \tab Hedgehog signaling pathway \cr
04350 \tab TGF-beta signaling pathway \cr
04360 \tab Axon guidance \cr
04370 \tab VEGF signaling pathway \cr
04510 \tab Focal adhesion \cr
04512 \tab ECM-receptor interaction \cr
04514 \tab Cell adhesion molecules (CAMs) \cr
04520 \tab Adherens junction \cr
04530 \tab Tight junction \cr
04540 \tab Gap junction \cr
04610 \tab Complement and coagulation cascades \cr
04612 \tab Antigen processing and presentation \cr
04614 \tab Renin-angiotensin system \cr
04620 \tab Toll-like receptor signaling pathway \cr
04630 \tab Jak-STAT signaling pathway \cr
04640 \tab Hematopoietic cell lineage \cr
04650 \tab Natural killer cell mediated cytotoxicity \cr
04660 \tab T cell receptor signaling pathway \cr
04662 \tab B cell receptor signaling pathway \cr
04664 \tab Fc epsilon RI signaling pathway \cr
04670 \tab Leukocyte transendothelial migration \cr
04710 \tab Circadian rhythm \cr
04720 \tab Long-term potentiation \cr
04730 \tab Long-term depression \cr
04740 \tab Olfactory transduction \cr
04742 \tab Taste transduction \cr
04810 \tab Regulation of actin cytoskeleton \cr
04910 \tab Insulin signaling pathway \cr
04912 \tab GnRH signaling pathway \cr
04914 \tab Progesterone-mediated oocyte maturation \cr
04916 \tab Melanogenesis \cr
04920 \tab Adipocytokine signaling pathway \cr
04930 \tab Type II diabetes mellitus \cr
04940 \tab Type I diabetes mellitus \cr
04950 \tab Maturity onset diabetes of the young \cr
05010 \tab Alzheimer's disease \cr
05012 \tab Parkinson's disease \cr
05014 \tab Amyotrophic lateral sclerosis (ALS) \cr
05110 \tab Vibrio cholerae infection \cr
05120 \tab Epithelial cell signaling in Helicobacter pylori infection \cr
05130 \tab Pathogenic Escherichia coli infection - EHEC \cr
05131 \tab Pathogenic Escherichia coli infection - EPEC \cr
05200 \tab Pathways in cancer \cr
05210 \tab Colorectal cancer \cr
05211 \tab Renal cell carcinoma \cr
05212 \tab Pancreatic cancer \cr
05213 \tab Endometrial cancer \cr
05214 \tab Glioma \cr
05215 \tab Prostate cancer \cr
05216 \tab Thyroid cancer \cr
05217 \tab Basal cell carcinoma \cr
05218 \tab Melanoma \cr
05219 \tab Bladder cancer \cr
05220 \tab Chronic myeloid leukemia \cr
05221 \tab Acute myeloid leukemia \cr
05222 \tab Small cell lung cancer \cr
05223 \tab Non-small cell lung cancer \cr
05310 \tab Asthma \cr
05320 \tab Autoimmune thyroid disease \cr
05322 \tab Systemic lupus erythematosus \cr
05330 \tab Allograft rejection \cr
05332 \tab Graft-versus-host disease \cr
05340 \tab Primary immunodeficiency 
}

}
\source{
http://www.genome.jp/kegg/
}
\references{
Kanehisa: Kyoto Encyclopaedia of Genes and Genomes
}
\examples{
\dontrun{
library(ddepn)
data(kegggraphs)
names(kegggraphs)
}
}
\keyword{datasets}
