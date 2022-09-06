"plot_mantel_1hierarchy" <- function(data12,
                                    data_node,
                                    mantel.node.size=10,
                                    mantel.label.size=1,
                                    edge.size.max.r_thre =0.6,
                                    edge.size.min.r_thre = 0.4,
                                    edge.size.max =10,
                                    edge.size.center= 5,
                                    edge.size.min = 1,
                                    edge.color.max.p_thre =0.05,
                                    edge.color.center.p_thre = 0.01,
                                    edge.color.min.p_thre = 0.001,
                                    edge.color.ultra ="grey90",
                                    edge.color.max ="red",
                                    edge.color.center="blue",
                                    edge.color.min = "darkgreen",
                                    note.text = c("Phenotype", "Mantel", "Phenotype"),
                                    note.text.color = c("#c51b7d", "#2f6661", "#c51b7d"),
                                    arrow.color.left ="#c51b7d",
                                    arrow.color.right = "#2f6661",
                                    label.cex=1,
                                    space.h=1.6,
                                    space.v=1.5,
                                    arrow.width = 0.5,
                                    arrow.size = 1.5,
                                    pbreaks = c(-Inf, 0.01, 0.05, Inf),
                                    plabels = c("< 0.01", "0.01 - 0.05", ">= 0.05"),
                                    rbreaks = c(-Inf,0, 0.2, 0.4, Inf),
                                    rlabels = c("< 0","0 - 0.2", "0.2 - 0.4", ">= 0.4"),
                                    size.r=c(0.5,1,2,4),
                                    color.p=random2color(5),
                                    node.label.color=NULL,
                                    vertex.receiver = seq(1,6),
                                    color.use=colorCustom(50,pal = "gygn")
) {
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

  net3 <- rbind(cbind(matrix(0, m1, m1), net2), ###原矩阵行+选择的行+添加的mantel行
                matrix(0,n1, m1 + n1))
  row.names(net3) <- c(row.names(net)[vertex.receiver],
                       row.names(net)[setdiff(1:m1,vertex.receiver)],
                       colnames(net2)[1:n1])
  colnames(net3) <- row.names(net3)


  color.use3 <- c(color.use[vertex.receiver],
                  color.use[setdiff(1:m1,vertex.receiver)],
                  rep("#ffffff", n1))

  color.use3.frame <- c(color.use[vertex.receiver],
                        color.use[setdiff(1:m1,vertex.receiver)],
                        color.use[seq(m1+1,m1+n1)])

  colnames(net3) <- row.names(net3)

  shape <- c(rep("circle", m), rep("circle", m1 - m),
             rep("circle", n1))

  g <- graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)

  edge.start <- ends(g, es = E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m, 1] <- 0
  coords[(m + 1):m1, 1] <- space.h
  coords[(m1 + 1):nrow(net3), 1] <- space.h/2
  coords[1:m, 2] <- seq(space.v, 0, by = -space.v/(m - 1))
  coords[(m + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 -
                                                            m - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1-1))

  coords_scale <- coords
  igraph::V(g)$size <- c(vertex.weight[seq(1,m)],vertex.weight[seq(m+1,m1)],
                         rep(mantel.node.size,n1))


  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]

  igraph::V(g)$label.cex <- label.cex
  igraph::V(g)$label.cex <-c(rep(label.cex,m1),rep(mantel.label.size,n1))

  if (is.null(node.label.color)) {
    igraph::V(g)$label.color <- color.use3.frame[igraph::V(g)]
  } else {
    igraph::V(g)$label.color <- node.label.color
  }

  E(g)$arrow.width <- arrow.width
  E(g)$arrow.size <- arrow.size


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

  label.dist <- c(rep(space.h * 2.8, m),
                  rep(space.h * 2.8, m1 - m),
                  rep(2.8, n1))
  label.locs <- c(rep(-pi, m),
                  rep(0, m1 - m),
                  rep(pi/2, n1-1),
                  rep(-pi/2, 1))

  text.pos <- cbind(c(-space.h/1.6, 0, space.h/1.6),
                    space.v - space.v/7)

  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle,
                           parameters = list(vertex.frame.color = 1,
                                             vertex.frame.width = 1))

  plot(g, edge.curved = 0, layout = coords_scale,
       margin = 0.2, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - m1)),
       vertex.label.degree = label.locs, vertex.label.dist = label.dist,
       vertex.label.family = "serif")

  text(text.pos, note.text, cex = 0.8,
       col = note.text.color,family="serif")
  arrow.pos1 <- c(-space.h/1.5, space.v - space.v/4, -0.04,
                  space.v - space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v - space.v/4, 0.04,
                  space.v - space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3],
                arrow.pos1[4], col = arrow.color.left, arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3],
                arrow.pos2[4], col = arrow.color.right, arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")

  title.pos = c(0, space.v)
  text(title.pos[1], title.pos[2], paste0("Micro-occurrence communication network"), cex = 1,family="serif")

}
