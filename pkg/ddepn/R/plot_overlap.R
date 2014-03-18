## add the edges from ddepn and color them
plot_overlap <- function(adjcore, adjX, overlap, symbols.plot, ig.layout, main, nameA, nameB, 
                         node.size1=25, node.size2=7, stimnodes=NULL, customtitle=NULL,
                         col_equal="black", col_new="blue", col_missing="red", vlabel.cex=1) {
  adjY <- adjcore
  ## default: all edges in adjY are colored gray
  mattoplot <- pmax(adjX,adjY)
  mattoplotsimple <- ifelse(mattoplot==0,0,1)
  igY <- graph.adjacency(mattoplotsimple)
  igedges <- gsub(" ","",as.matrix(print.igraph.es(E(igY))))
  ecols <- rep("gray",length(igedges))
  ewidths <- rep(1,length(igedges))
  eltys <- rep("solid",length(igedges))
  ## traverse all nodes that are in both networks
  for(ov1 in overlap) {
    for(ov2 in overlap) {
      ## set edge type
      newtype <- adjX[ov1,ov2]
      oldtype <- adjY[ov1,ov2]
      ee <- paste(ov1,ov2,sep="->")
      wee <- which(igedges==ee)
      # no edge in new net
      if(newtype==0) {
        coltmp <- col_missing
        ltytmp <- "solid"
        if(oldtype==0) { ## true negative
          coltmp <- "white" # don't show an edge, none there
          next
        }
        if(oldtype==2) { ## inhibition lost
          ltytmp <- "dashed"
        }  			
      } else if(newtype==1) {
        coltmp <- col_new
        ltytmp <- "solid"
        if(oldtype==1) { ## edge is the same and correct inferred
          coltmp <- col_equal
        }
        if(oldtype==2) { ## type switched
          coltmp <- "orange"
        }				
      } else if(newtype==2) {
        coltmp <- col_new #"blue"
        ltytmp <- "dashed"
        if(oldtype==1) { ## type switched
          coltmp <- "orange"
        }
        if(oldtype==2) { ## edge is the same and correct inferred
          coltmp <- col_equal
        }				
      }
      print(paste(coltmp,ltytmp,3))
      ecols[wee] <- coltmp # change the color
      eltys[wee] <- ltytmp # change the line type
      ewidths[wee] <- 3 # set as thick line, to mark that the edge was included in the net inference
      
    }
  }
  ## make reciprocated edges be curved
  ecurved <- rep(0,length(igedges))
  ed <- which(mattoplot!=0,arr.ind=TRUE)
  erecipn <- NULL
  for(ei in 1:nrow(ed)) {
    edd <- ed[ei,]
    eddg <- mattoplot[edd[1],edd[2]]
    eddgr <- mattoplot[edd[2],edd[1]]
    if(eddgr!=0)
      erecipn <- c(erecipn,paste(rownames(mattoplot)[edd[1]],colnames(mattoplot)[edd[2]],sep="->"))
  }
  ecurved[match(erecipn,igedges)] <- 0.3
  E(igY)$curved <- ecurved
  
  ig.nodes <- gsub(" ","",as.matrix(print.igraph.vs(V(igY))))
  print("Setting graph attributes...")
  if(is.null(stimnodes)) {
    V(igY)$color <- "gray"
  } else {
    vcols <- rep("gray",length(ig.nodes))
    vcols[match(stimnodes, ig.nodes)] <- "red"
    V(igY)$color <- vcols
  }
  V(igY)$shape <- rep("crectangle",length(ig.nodes))
  V(igY)$size <- rep(node.size1,length(ig.nodes))
  V(igY)$size2 <- rep(node.size2,length(ig.nodes))
  V(igY)$label.cex <- vlabel.cex
  E(igY)$color <- ecols
  E(igY)$lty <- eltys
  E(igY)$width <- ewidths
  
  if(!is.null(customtitle)) {
    main <- customtitle
  } else {
    main <- c(paste(nameA, "vs.", nameB),main)
  }
  plot(igY,vertex.label=symbols.plot,layout=ig.layout,main=main)
  cls <- c("gray",col_equal,"blue","red","orange","black")
  lts <- c("solid","solid","solid","solid","solid","dashed")
  lwd <- c(1,3,3,3,3)
  legs <- c("only in prior","equal","added","deleted","switched type","inhibition")
  legend("bottomright", xpd=NA, legend=legs, col=cls, lty=lts, lwd=lwd)
}