"cor_mantel" <- function(edmat,ctmat,
                         method="pearson",
                         permutations=9999,
                         filter=TRUE,ed.p=0.05,ed.r=0.4,
                         mantel_select=list(Spec01 = 1:7)) {

  occor.r<-WGCNA::corAndPvalue(ctmat,method = 'pearson')$cor
  occor.p<-WGCNA::corAndPvalue(ctmat,method = 'pearson')$p
  diag(occor.r) <- 0

  mtadj<-multtest::mt.rawp2adjp(unlist(occor.p),proc='Bonferroni')
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  adjp<-adpcor
  adjp<-matrix(adjp,dim(occor.p)[2])
  row.names(adjp)<-row.names(occor.p)
  colnames(adjp)<-colnames(occor.p)
  occor.p<-adjp

  if (filter) occor.r[occor.p>ed.p|abs(occor.r)<ed.r] = 0

  mantel <- mantel_calc( edmat, ctmat,
                         method=method,
                         permutations=permutations,
                         select_col =mantel_select)

  return(list(occor.r=occor.r,mantel=mantel,occor.p=occor.p))
}


