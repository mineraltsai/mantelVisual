"plot_mantel_2hierarchy" <- function(data_test,
                                     data_node, ###计算节点大小
                                     margin=0.2,
                                     mantel.node.size=10,
                                     mantel.label.size=0.5,
                                     label.cex=1,
                                     edge.size.max.r_thre = 0.6,
                                     edge.size.min.r_thre = 0.4,
                                     edge.size.max = 8,
                                     edge.size.center= 4,
                                     edge.size.min = 2,
                                     note.text = c("Phenotype", "Mantel", "Phenotype"),
                                     note.text.color = c("#c51b7d", "#2f6661", "#c51b7d"),
                                     arrow.color.left ="#c51b7d",
                                     arrow.color.right = "#2f6661",
                                     space.h=1.6,
                                     space.v=1.5,
                                     arrow.width = 0,
                                     arrow.size = 0,
                                     pbreaks = c(-Inf,0.001, 0.01, 0.05, Inf),
                                     plabels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05"),
                                     color.p=c("grey90","red","blue","yellow"),
                                     rbreaks = c(-Inf,0, 0.4, 0.8, Inf),
                                     rlabels = c("< 0","0 - 0.2", "0.2 - 0.8", ">= 0.8"),
                                     size.r=c(0.5,2,5,10),
                                     node.label.color=NULL,
                                     vertex.receiver = seq(1,6),
                                     color.use=colorCustom(50,pal = "gygn"))
{
  library(tidyverse)

  corr<-data12$occor.r
  d2<-data12$mantel
  corp<-data12$occor.p

  d3<-subset(d2,select=c(1:3))
  d4<-subset(d2,select=c(1:2,4))

  d5<-tidyr::spread(d3,key = "spec",value = r,drop=FALSE)
  d5<-tibble::column_to_rownames(d5,var = "env")

  d6<-tidyr::spread(d4,key = "spec",value = p,drop=FALSE)
  d6<-tibble::column_to_rownames(d6,var = "env")


  mantelr<-d5%>%as.matrix()%>%abs()
  mantelp<-d6%>%as.matrix()%>%abs()

  netc<-matrMerge(corr,mantelr)
  netp<-matpMerge(corp,mantelp)
  net<-corr%>%abs()
  rownames(net)<-colnames(net)

  library(micro2eco)
  data<-data_node
  sum<-list()
  mean<-list()
  colnum<-length(colnames(data))
  for (k in 1:colnum) {
    data[,k]<-micro2eco::normalize(data[,k])$x
    sum[k]<-sum(data[,k])
    mean[k]<-mean(data[,k])
  }

  vertex.weight<-sum%>%as.numeric()
  vertex.weight.max <- max(vertex.weight)

  if (length(unique(vertex.weight)) == 1) {
    vertex.size.max <- 5
  }else {
    vertex.size.max <- vertex.weight%>%max()
    vertex.weight.max<- vertex.size.max
  }

  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max

  net2 <- netc
  m <- length(vertex.receiver)
  m1 <- nrow(net2)
  n1 <- ncol(net2)

  net3 <- rbind(cbind(matrix(0, m1, m1),net2), ###原矩阵行+选择的行+添加的mantel行
                matrix(0, n1, m1 + n1))


  row.names(net3) <- c(row.names(corr),
                       colnames(net2)[1:n1])
  colnames(net3) <- row.names(net3)


  color.use3 <- c(
    color.use[1:nrow(corr)],
    rep("#ffffff", n1))

  color.use3.frame <- c(
    color.use[1:nrow(corr)],
    color.use[seq(m1+1,m1+n1)])

  colnames(net3) <- row.names(net3)

  shape <- c(
    rep("circle", m1),
    rep("circle", n1))

  g <- graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)

  edge.start <- ends(g, es = E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m1, 1] <- 0
  coords[(m1 + 1):nrow(net3), 1] <- space.h

  coords[1:m1, 2] <- seq(space.v, 0, by = -space.v/(m1 - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v*9/10, space.v*1/10, length.out=n1)

  coords_scale <- coords
  igraph::V(g)$size <- c(vertex.weight[seq(1,m1)],
                         rep(5,n1))


  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]


  igraph::V(g)$label.cex <-c(rep(1,m1),
                             rep(1,n1))

  if (is.null(node.label.color)) {
    igraph::V(g)$label.color <- color.use3.frame[igraph::V(g)]
  } else {
    igraph::V(g)$label.color <- node.label.color
  }

  E(g)$arrow.width <- 0
  E(g)$arrow.size <- 0



  netp<-t(netp)

  rd = cut(netc, breaks = rbreaks ,
           labels = rlabels)

  r_num<-length(levels(rd))
  size.use2<-size.r
  names(size.use2) <-levels(rd)
  size.use2<-size.use2%>%data.frame()
  size.use2<-rownames_to_column(size.use2,var = "color")

  rd<-data.frame(rd = cut(netc, breaks = rbreaks,
                          labels = rlabels))

  rd$rd<-as.character(rd$rd)
  rd$color<-rd$rd
  for (k in 1:r_num) {
    rd$color[which(match(rd$rd,size.use2$color[k])==1)]<-size.use2$.[k]
  }
  rd<-rd$color
  rrd<-matrix(rd,nrow(netc),ncol(netc))
  colnames(rrd)<-colnames(netc)
  rownames(rrd)<-rownames(netc)
  E(g)$weight<-rrd%>%t()
  E(g)$width <- E(g)$weight





  pd = cut(netp, breaks = pbreaks,
           labels = plabels)

  p_num<-length(levels(pd))
  color.use2<-color.p
  names(color.use2) <-levels(pd)
  color.use2<-color.use2%>%data.frame()
  color.use2<-rownames_to_column(color.use2,var = "color")

  pd<-data.frame(pd = cut(netp, breaks = pbreaks,
                          labels = plabels))
  pd$pd<-as.character(pd$pd)
  pd$color<-pd$pd
  for (k in 1:p_num) {
    pd$color[which(match(pd$pd,color.use2$color[k])==1)]<-color.use2$.[k]
  }
  pd<-pd$color
  ppd<-matrix(pd,nrow(netp),ncol(netp))
  colnames(ppd)<-colnames(netp)
  rownames(ppd)<-rownames(netp)

  E(g)$color<-ppd



  label.dist<-NULL
  label.locs<-NULL

  label.dist <- c(rep(space.h * 2.8, m1),
                  rep(2.8, n1))
  label.locs <- c(
    rep(-pi, m1),
    rep(0, n1-1),
    rep(0, 1))

  text.pos <- cbind(c(-space.h/1.6,  space.h/1.6),
                    space.v - space.v/7)

  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle,
                           parameters = list(vertex.frame.color = 1,
                                             vertex.frame.width = 1))

  plot(g, edge.curved = 0, layout = coords_scale,
       margin = margin, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - m1)),
       vertex.label.degree = label.locs, vertex.label.dist = label.dist,
       vertex.label.family = "serif")

  text(text.pos, note.text, cex = 0.8,
       col = note.text.color,family="serif")
  arrow.pos1 <- c(-space.h/1.5, space.v - space.v/4, space.h/1.5, space.v - space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3],
                arrow.pos1[4], col = arrow.color.left, arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")

  title.pos = c(0, space.v)
  text(title.pos[1], title.pos[2], paste0("Micro-occurrence network"), cex = 1,family="serif")



}
