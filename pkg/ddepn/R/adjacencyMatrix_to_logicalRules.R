# function to convert adjacency matrices of 'ddepn' consensus nets to logical
# activation/inhibition rules of network component, i.e. target, wiring for 
# further analysis with the R package 'BoolNet', e.g. for perturbation simulations,
# according to its 'loadNetwork'-function format
# 
# Author: Silvia von der Heyde
###################################################################################

adjacencyMatrix_to_logicalRules<-function(adjMatrix, outfile){

## get targets from input adjacency matrix
targets<-colnames(adjMatrix)  # should equal rownames
if(any(targets!=rownames(adjMatrix))){
  stop("Error in adjacency matrix: rownames don't equal colnames!")
}

## initialise net description for 'BoolNet' 
sink(outfile,append=TRUE)
cat("targets , factors\n")
sink()

n<-length(targets)
for(i in 1:n){
    tar<-targets[i]
    colv<-adjMatrix[,i]
    w.1<-which(colv==1) # activating factors (at least one necessary for activation)
    w.2<-which(colv==2) # inhibiting factors (at least one necessary for inhibition)
    w.0<-which(colv==0)
    if(length(w.0)==n){ # no influence from other factors --> constant
         sink(outfile,append=TRUE)
         cat(paste(tar,",",sep=""),tar,"\n")
         sink()
    }
    if(length(w.1)==0 & length(w.2)>0){ # no activators
         tar.inh<-rownames(adjMatrix)[w.2]    
         sink(outfile,append=TRUE)
         #  tar & !(i1 | i2 | ...)
         cat(paste(tar,",",sep=""),tar,"&",paste("!(",paste(tar.inh,collapse="|"),")",sep=""),"\n")
         sink()
    }
    if(length(w.1)>0 & length(w.2)==0){ # no inhibitors
         tar.act<-rownames(adjMatrix)[w.1]    
         sink(outfile,append=TRUE)
         # tar | a1 | a2 | ...
         cat(paste(tar,",",sep=""),paste(c(tar,tar.act),collapse="|"),"\n")
         sink()
    }
    if(length(w.1)>0 & length(w.2)>0){ # inhibitors AND activators
         tar.act<-rownames(adjMatrix)[w.1]    
         tar.inh<-rownames(adjMatrix)[w.2]  
         sink(outfile,append=TRUE)
         # (tar | a1 | a2 | ...) & !(i1 | i2 | ...)
         cat(paste(tar,",",sep=""),paste("(",paste(c(tar,tar.act),collapse="|"),")",sep=""),"&",paste("!(",paste(tar.inh,collapse="|"),")",sep=""),"\n")
         sink()
    }
}

}