crossover <- function(phi1, phi2) {
  crossnode <- sample(1:nrow(phi1),1)
  phi1.tmp <- phi1
  # change line
  phi1[crossnode,] <- phi2[crossnode,]
  phi2[crossnode,] <- phi1.tmp[crossnode,]
  # now change the column except for the position phiX[crossnode,crossnode], since it is already changed
  phi1[-crossnode,crossnode] <- phi2[-crossnode,crossnode]
  phi2[-crossnode,crossnode] <- phi1.tmp[-crossnode,crossnode]
  return(list(phi1,phi2))
}
