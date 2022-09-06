
"mantel_calc" <- function(spec,param,
                          method="pearson",
                          permutations=9999,
                          select_col=list( firm = 2:2,cya= 1:1,

                                           pro=3:3)) {
  options(warn = -1)
  cat("\n","\n",'method: "pearson", "spearman" or "kendall"')
  datad<-data.frame()
  mantel.result<-data.frame()
  for (k in 1:length(select_col)) {
    select_col1<-select_col[[k]]
    spec.name<-names(select_col)[[k]]
    spec[,select_col1][which(spec[,select_col1]==0)]<-1
    veg.dist <- vegan::vegdist(spec[,select_col1]) # Bray-Curtis

    for (t in 1:ncol(param)) {
      env.dist <- vegan::vegdist(scale(param[,t]), "euclid")
      mantel.res<-vegan::mantel(veg.dist, env.dist,
                                permutations=permutations,
                                method=method)
      mantel.res.r<-mantel.res$statistic
      mantel.res.p<-mantel.res$signif
      dad<-data.frame(spec=spec.name,
                      env=colnames(param)[t],
                      r=mantel.res.r,
                      p=mantel.res.p)

      {
        datad<-rbind(datad,dad)
      }

    }

  }
  datad$spec<-factor(datad$spec,levels = names(select_col))

  return(datad)
}


