"plot_mantel_heatmap" <- function(data12,
                           edge.curve.min = -0.3,
                           edge.curve.max = 0.3,add.sig=TRUE,
                           layout = c("lower","upper"),
                           grid.background =TRUE,
                           mantel.node.pos.para=TRUE, ###false,不平行，则需要调整下面参数
                           mantel.node.x.pos=c(5,7,9,11),
                           mantel.node.y.pos=c(13,11,9,7),
                           adjust.diag=TRUE, ###是否调整对角线距离
                           adjust.diag.para=TRUE, ###是否平行调整对角线距离
                           adjust.diag.para.dist=1,  ###是否平行调整对角线距离为1
                           adjust.node = c(5,7),
                           adjust.diag.dist = c(3,4),
                           diag.node.size = 5,
                           heatmap.node.size.ratio=15,
                           mantel.node.size.para=TRUE,
                           mantel.node.size=5,
                           mantel.node.size.custom=rep(5,4),  ###5为大小，4为数量
                           heatmap.shape = 1,
                           diag.node.shape.order = 1,
                           mantel.node.shape.order = c(1,2,3,4),
                           adjust.text.x = 0.8,
                           adjust.text.y = 0.8,
                           grid.color = "grey90",
                           diag.frame.color = "red",
                           mantel.frame.color="purple",
                           space.h=1.6,
                           space.v=1.5,
                           node.color.postive = "blue",
                           node.color.mid = "white",
                           node.color.negative = "red",
                           diag.node.color = "red",
                           mantel.node.color = "purple",
                           label.cex=1,
                           edge.size.max = 8,
                           edge.size.center= 4,
                           edge.size.min = 2,
                           arrow.width = 0,
                           arrow.size = 0,
                           pbreaks = c(-Inf,0.001, 0.01, 0.05, Inf),
                           plabels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05"),
                           color.p=c("grey90","red","blue","yellow"),
                           rbreaks = c(-Inf,0, 0.4, 0.8, Inf),
                           rlabels = c("< 0","0 - 0.2", "0.2 - 0.8", ">= 0.8"),
                           size.r=c(0.5,2,5,10),
                           node.label.color=NULL) {

  library(tidyverse)
  if (grid.background) {
    "mantel_data" <- function(data12) {

    corr<-data12$occor.r
    corp<-data12$occor.p
    mantel<-data12$mantel

    d3<-subset(mantel,select=c(1:3))
    d4<-subset(mantel,select=c(1:2,4))

    d5<-tidyr::spread(d3,key = "spec",value = r,drop=FALSE)
    d5<-tibble::column_to_rownames(d5,var = "env")

    d6<-tidyr::spread(d4,key = "spec",value = p,drop=FALSE)
    d6<-tibble::column_to_rownames(d6,var = "env")
    mantelr<-d5%>%as.matrix()%>%abs()
    mantelp<-d6%>%as.matrix()%>%abs()

    netc<-matrMerge(corr,mantelr)
    netp<-matpMerge(corp,mantelp)


    return(list(corr=corr,corp=corp,
                diagcorr=netc,diagp=netp))
  }

  mantel_data<-mantel_data(data12)
  net3<-mantel_data$corr%>%abs()
  corr<-mantel_data$corr
  corp<-mantel_data$corp

  "coord_matrix" <- function(net3)
  {
    mat1<-net3
    diag1<-length(row.names(mat1))
    diag<-diag1+1
    for (kk in 2:(nrow(net3)-1)) {
      matk<-net3[kk:nrow(net3),kk:nrow(net3)]

      {
        mat1<-rbind(cbind(mat1,matrix(0, ncol(mat1), nrow(matk))),
                    cbind(matrix(0, nrow(matk), ncol(mat1)),matk))
        diag_num<-length(row.names(mat1))+1
        diag<-rbind(diag,diag_num)
      }
    }
    matend<-net3[nrow(net3),nrow(net3)]
    mat1<-rbind(cbind(mat1,matrix(0, ncol(mat1), 1)),
                cbind(matrix(0, 1, ncol(mat1)),matend))
    return(list(mat1=mat1,diag=diag))
  }

  mat11<-coord_matrix(net3)
  mat1<-mat11$mat1
  mat1[mat1 != 0]<-0
  diag<-mat11$diag
  diag<-rbind(1,diag)

  ###添加对角线相关性
  netc<-mantel_data$diagcorr
  netp<-mantel_data$diagp

  mat2<-rbind(cbind(mat1,matrix(0, nrow(mat1), ncol(netc))),
              matrix(0, ncol(netc),  nrow(mat1)+ncol(netc)))

  for (tt in 1:length(diag)) {
    mat2[diag[tt],((nrow(mat1)+1):(nrow(mat1)+ncol(netc)))]<-netc[tt,]
  }


  row.names(mat2)[seq(nrow(mat1)+1,nrow(mat2))]<-colnames(netc)
  colnames(mat2)[seq(nrow(mat1)+1,nrow(mat2))]<-colnames(netc)

  diag_l<-diag%>%as.numeric()

  colnames(mat2)[diag_l]<-row.names(corr)

  pp<-corrposcalcu(corr=corr,p.mat=corp,method = "square",type = "lower")


  g <- graph_from_adjacency_matrix(mat2, mode = "undirected",
                                   weighted = T)

  coord_layout <- pp$corrPos[,3:4]

  xlim<-coord_layout$x%>%min()
  xmax<-coord_layout$x%>%max()
  ylim<-coord_layout$y%>%min()
  ymax<-coord_layout$y%>%max()

  add.node<-ncol(netc)


  if (mantel.node.pos.para) {
    xscale<-seq(as.integer(nrow(netc)/2),nrow(netc),length.out=add.node)
    yscale<- ymax+7-xscale
  } else {
    xscale<- mantel.node.x.pos
    yscale<- mantel.node.y.pos
  }

  add.coord<-data.frame(x=xscale,
                        y=yscale)

  data_coord<-rbind(coord_layout,add.coord)

  sub_net_layout <- data_coord%>%as.matrix()

  if (adjust.diag) {
    if (adjust.diag.para){
      sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1] <- sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1]+adjust.diag.para.dist
    } else {
      sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1] <- sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1]+adjust.diag.para.dist

      adjust.data <- sub_net_layout[c(diag), 1]%>%data.frame()

      diag_data<-data.frame(row=c(diag),
                            x=adjust.data$.)

      dd<-diag_data[adjust.node,]
      dd$x<-dd$x+adjust.diag.dist
      sub_net_layout[dd$row, 1] <-dd$x

    }
  } else {
    sub_net_layout<-sub_net_layout
  }

  if (mantel.node.size.para) {
    igraph::V(g)$weight <- c(pp$corrPos$corr%>%abs()%>%as.matrix()*heatmap.node.size.ratio,
                             rep(mantel.node.size,add.node))
  } else {
    igraph::V(g)$weight <- c(pp$corrPos$corr%>%abs()%>%as.matrix()*heatmap.node.size.ratio,
                             mantel.node.size.custom)
  }



  V(g)$weight[diag] <-diag.node.size

  igraph::V(g)$size <- V(g)$weight

  ###节点形状
  cat("\n","shape: ",shapes(),sep = "\n")
  diag.node.shape = shapes()[diag.node.shape.order]
  mantel.node.shape = shapes()[mantel.node.shape.order]
  ##定义所有节点形状
  shape_sel<-data.frame(shape=shapes()[heatmap.shape],
                        num=seq(1,nrow(mat2)))
  ##定义对角线节点形状
  shape_sel$shape[which(shape_sel$num %in% diag)] <-diag.node.shape
  ##定义mantel节点形状
  if (is.null(mantel.node.shape.order)) {
    shape_sel$shape[(nrow(mat2)-add.node+1):nrow(mat2)] <-shapes()[1]
  } else {
    shape_sel$shape[(nrow(mat2)-add.node+1):nrow(mat2)] <-mantel.node.shape

  }
  shapes<-shape_sel$shape%>%as.character()



  label.dist<-NULL
  label.locs<-NULL

  space.h<-adjust.text.x
  space.v<-adjust.text.y

  diag.label.dist<-space.h * 2.8
  mantel.label.dist<-space.h * 2.8
  label.dist<-rep(0,nrow(mat2))
  label.dist[2:nrow(corr)]<-space.h * 2.8
  label.dist[diag_l]<-diag.label.dist
  label.dist[seq(nrow(mat2)-add.node+1,nrow(mat2))]<-mantel.label.dist*(-1)


  label.locs<-rep(0,nrow(mat2))
  label.locs[2:nrow(corr)]<-pi*(-1)
  label.locs[diag_l]<-pi*(-1)
  label.locs[seq(nrow(mat2)-add.node+1,nrow(mat2))]<-pi

  ###边的颜色
  {netp<-t(netp)
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

    E(g)$color<-ppd}
  ###节点的颜色
  {node.display<- pp$corrPos[,5:6]

    color1_num<-node.display[which(node.display$corr<0),]%>%rownames()%>%length()
    color3_num<-node.display[which(node.display$corr>0),]%>%rownames()%>%length()

    node.display$color<-1
    node.display$frame.color<-1
    node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$color<-colorRampPalette(c(node.color.postive,node.color.mid),bias=1)(color1_num)
    node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$color<-colorRampPalette(c(node.color.mid,node.color.negative),bias=1)(color3_num)

    node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$frame.color<-node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$color
    node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$frame.color<-node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$color

    node.display[which(node.display$corr==0),]$color<-diag.node.color

    add.node.color<-rep(mantel.node.color,add.node)
    add.node.data<-data.frame(corr=rep(1,add.node),
                              p.value=rep(0.05,add.node),
                              color=add.node.color,
                              frame.color=add.node.color)


    node.display_comb<-rbind(node.display,add.node.data)

    igraph::V(g)$color <- node.display_comb$color
    igraph::V(g)$frame.color <- node.display_comb$frame.color
  }
  ###grid节点大小
  V(g)$frame.size<-max(V(g)$size)

  newsize<-V(g)$frame.size
  newsize[diag] <-V(g)$size[diag]
  newsize[(nrow(mat2)-add.node+1):nrow(mat2)] <-0

  newcolor<-rep(grid.color,nrow(mat2))
  newcolor[diag] <-diag.frame.color
  newcolor[(nrow(mat2)-add.node+1):nrow(mat2)] <-mantel.frame.color
  #??Plot.igraph
  ###add background grid
  g1<-g
  g2<-g1%>%add_vertices(nrow(mat2),
                        frame.color=newcolor,
                        name="",
                        frame.size=newsize,
                        size=newsize,
                        weight=newsize)


  sub_net_layout2<-rbind(sub_net_layout,sub_net_layout)
  layout=layout

  if (add.sig) {
    corp<-pp$corrPos[,1:6]
    corp$p.value<-corp$p.value%>%as.double()
    corp[which(corp$p.value<=0.001 & corp$p.value >0),]$p.value <-"***"
    corp[which(corp$p.value<=0.01 & corp$p.value >0.001),]$p.value <-"**"
    corp[which(corp$p.value<0.05 & corp$p.value >0.01),]$p.value <-"*"
    corp[which(corp$p.value>=0.05 | corp$p.value==0),]$p.value<-""
    corp<-corp$p.value%>%as.matrix()

    newp<-corp
    newp[diag] <-""
    newp[(nrow(mat2)-add.node+1):nrow(mat2)] <-""
    vertex_attr(g2)$name[seq(nrow(mat2)+1,nrow(mat2)+nrow(mat2))]<-newp


    label.dist1<-label.dist
    label.dist1[seq(1,nrow(mat2))]<-0
    label.dist1[diag] <- 0
    label.dist1[(nrow(mat2)-add.node+1):nrow(mat2)] <-0
    label.dist2<-c(label.dist,label.dist1)
  } else {
    label.dist2<-label.dist
  }

  if (layout == "lower") {

    edge.curved1<- matrix(rep(1,add.node),ncol(netp),nrow(netp))
    for (jj in 1:add.node) {
      order <- jj*2
      edge.curved<-rep(1,nrow(corr))
      edge.curved[order]<-0
      edge.curved[1:order]<-seq(edge.curve.max,0,length.out=order)
      edge.curved[order:length(edge.curved)]<-seq(0,edge.curve.min,length.out=length(edge.curved)-order+1)
      edge.curved1[,jj]<-edge.curved
    }

    edge.curved1<-t(edge.curved1)
    E(g1)$curved<-edge.curved1

    label.locs1<-label.locs*(-1)
    label.dist1<-label.dist*(-1)

    sub_net_layout1<-sub_net_layout2%>%data.frame()
    sub_net_layout1$x<-nrow(corr)+3-sub_net_layout1$x
    sub_net_layout1$y<-nrow(corr)+3-sub_net_layout1$y
    sub_net_layout1<-sub_net_layout1%>%as.matrix()

    label.dist1<-c(label.dist1,label.dist2[(nrow(mat2)+1):(nrow(mat2)+nrow(mat2))])
    ##vertex.attributes(g2)
    plot(g2,edge.curved = edge.curved1,
         vertex.shape=shapes,
         vertex.label=vertex_attr(g2)$name,
         layout=sub_net_layout1,
         vertex.label.degree = label.locs1,
         vertex.label.dist = label.dist1,
         vertex.label.family = "serif")
  }

  if (layout == "upper") {##vertex.attributes(g2)

    edge.curved1<- matrix(rep(1,add.node),ncol(netp),nrow(netp))
    for (jj in 1:add.node) {
      order <- jj*2
      edge.curved<-rep(1,nrow(corr))
      edge.curved[order]<-0
      edge.curved[1:order]<-seq(edge.curve.max,0,length.out=order)
      edge.curved[order:length(edge.curved)]<-seq(0,edge.curve.min,length.out=length(edge.curved)-order+1)
      edge.curved1[,jj]<-edge.curved
    }

    edge.curved1<-t(edge.curved1)
    E(g1)$curved<-edge.curved1

    plot(g2,edge.curved = edge.curved1,
         vertex.shape=shapes,
         vertex.label=vertex_attr(g2)$name,
         layout=sub_net_layout2,
         vertex.label.degree = label.locs,
         vertex.label.dist = label.dist2,
         vertex.label.family = "serif")}} else {

           "mantel_data" <- function(data12) {

           corr<-data12$occor.r
           corp<-data12$occor.p
           mantel<-data12$mantel

           d3<-subset(mantel,select=c(1:3))
           d4<-subset(mantel,select=c(1:2,4))

           d5<-tidyr::spread(d3,key = "spec",value = r,drop=FALSE)
           d5<-tibble::column_to_rownames(d5,var = "env")

           d6<-tidyr::spread(d4,key = "spec",value = p,drop=FALSE)
           d6<-tibble::column_to_rownames(d6,var = "env")
           mantelr<-d5%>%as.matrix()%>%abs()
           mantelp<-d6%>%as.matrix()%>%abs()

           netc<-matrMerge(corr,mantelr)
           netp<-matpMerge(corp,mantelp)


           return(list(corr=corr,corp=corp,
                       diagcorr=netc,diagp=netp))
         }

         mantel_data<-mantel_data(data12)
         net3<-mantel_data$corr%>%abs()
         corr<-mantel_data$corr
         corp<-mantel_data$corp

         "coord_matrix" <- function(net3)
         {
           mat1<-net3
           diag1<-length(row.names(mat1))
           diag<-diag1+1
           for (kk in 2:(nrow(net3)-1)) {
             matk<-net3[kk:nrow(net3),kk:nrow(net3)]

             {
               mat1<-rbind(cbind(mat1,matrix(0, ncol(mat1), nrow(matk))),
                           cbind(matrix(0, nrow(matk), ncol(mat1)),matk))
               diag_num<-length(row.names(mat1))+1
               diag<-rbind(diag,diag_num)
             }
           }
           matend<-net3[nrow(net3),nrow(net3)]
           mat1<-rbind(cbind(mat1,matrix(0, ncol(mat1), 1)),
                       cbind(matrix(0, 1, ncol(mat1)),matend))
           return(list(mat1=mat1,diag=diag))
         }

         mat11<-coord_matrix(net3)
         mat1<-mat11$mat1
         mat1[mat1 != 0]<-0
         diag<-mat11$diag
         diag<-rbind(1,diag)

         ###添加对角线相关性
         netc<-mantel_data$diagcorr
         netp<-mantel_data$diagp

         mat2<-rbind(cbind(mat1,matrix(0, nrow(mat1), ncol(netc))),
                     matrix(0, ncol(netc),  nrow(mat1)+ncol(netc)))

         for (tt in 1:length(diag)) {
           mat2[diag[tt],((nrow(mat1)+1):(nrow(mat1)+ncol(netc)))]<-netc[tt,]
         }


         row.names(mat2)[seq(nrow(mat1)+1,nrow(mat2))]<-colnames(netc)
         colnames(mat2)[seq(nrow(mat1)+1,nrow(mat2))]<-colnames(netc)

         diag_l<-diag%>%as.numeric()

         colnames(mat2)[diag_l]<-row.names(corr)

         pp<-corrposcalcu(corr=corr,p.mat=corp,method = "square",type = "lower")


         g <- graph_from_adjacency_matrix(mat2, mode = "undirected",
                                          weighted = T)

         coord_layout <- pp$corrPos[,3:4]

         xlim<-coord_layout$x%>%min()
         xmax<-coord_layout$x%>%max()
         ylim<-coord_layout$y%>%min()
         ymax<-coord_layout$y%>%max()

         add.node<-ncol(netc)

         if (mantel.node.pos.para) {
           xscale<-seq(as.integer(nrow(netc)/2),nrow(netc),length.out=add.node)
           yscale<- ymax+7-xscale
         } else {
           xscale<- mantel.node.x.pos
           yscale<- mantel.node.y.pos
         }

         add.coord<-data.frame(x=xscale,
                               y=yscale)

         data_coord<-rbind(coord_layout,add.coord)

         sub_net_layout <- data_coord%>%as.matrix()

         if (adjust.diag) {
           if (adjust.diag.para){
             sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1] <- sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1]+adjust.diag.para.dist
           } else {
             sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1] <- sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1]+adjust.diag.para.dist

             adjust.data <- sub_net_layout[c(diag), 1]%>%data.frame()

             diag_data<-data.frame(row=c(diag),
                                   x=adjust.data$.)

             dd<-diag_data[adjust.node,]
             dd$x<-dd$x+adjust.diag.dist
             sub_net_layout[dd$row, 1] <-dd$x

           }
         } else {
           sub_net_layout<-sub_net_layout
         }

         if (mantel.node.size.para) {
           igraph::V(g)$weight <- c(pp$corrPos$corr%>%abs()%>%as.matrix()*heatmap.node.size.ratio,
                                    rep(mantel.node.size,add.node))
         } else {
           igraph::V(g)$weight <- c(pp$corrPos$corr%>%abs()%>%as.matrix()*heatmap.node.size.ratio,
                                    mantel.node.size.custom)
         }



         V(g)$weight[diag] <-diag.node.size

         igraph::V(g)$size <- V(g)$weight

         ###节点形状
         cat("\n","shape: ",shapes(),sep = "\n")
         diag.node.shape = shapes()[diag.node.shape.order]
         mantel.node.shape = shapes()[mantel.node.shape.order]
         ##定义所有节点形状
         shape_sel<-data.frame(shape=shapes()[heatmap.shape],
                               num=seq(1,nrow(mat2)))
         ##定义对角线节点形状
         shape_sel$shape[which(shape_sel$num %in% diag)] <-diag.node.shape
         ##定义mantel节点形状
         if (is.null(mantel.node.shape.order)) {
           shape_sel$shape[(nrow(mat2)-add.node+1):nrow(mat2)] <-shapes()[1]
         } else {
           shape_sel$shape[(nrow(mat2)-add.node+1):nrow(mat2)] <-mantel.node.shape

         }
         shapes<-shape_sel$shape%>%as.character()



         label.dist<-NULL
         label.locs<-NULL

         space.h<-adjust.text.x
         space.v<-adjust.text.y

         diag.label.dist<-space.h * 2.8
         mantel.label.dist<-space.h * 2.8
         label.dist<-rep(0,nrow(mat2))
         label.dist[2:nrow(corr)]<-space.h * 2.8
         label.dist[diag_l]<-diag.label.dist
         label.dist[seq(nrow(mat2)-add.node+1,nrow(mat2))]<-mantel.label.dist*(-1)


         label.locs<-rep(0,nrow(mat2))
         label.locs[2:nrow(corr)]<-pi*(-1)
         label.locs[diag_l]<-pi*(-1)
         label.locs[seq(nrow(mat2)-add.node+1,nrow(mat2))]<-pi

         ###边的颜色
         {netp<-t(netp)
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

           E(g)$color<-ppd}
         ###节点的颜色
         {node.display<- pp$corrPos[,5:6]

           color1_num<-node.display[which(node.display$corr<0),]%>%rownames()%>%length()
           color3_num<-node.display[which(node.display$corr>0),]%>%rownames()%>%length()

           node.display$color<-1
           node.display$frame.color<-1
           node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$color<-colorRampPalette(c(node.color.postive,node.color.mid),bias=1)(color1_num)
           node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$color<-colorRampPalette(c(node.color.mid,node.color.negative),bias=1)(color3_num)

           node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$frame.color<-node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$color
           node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$frame.color<-node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$color

           node.display[which(node.display$corr==0),]$color<-diag.node.color

           add.node.color<-rep(mantel.node.color,add.node)
           add.node.data<-data.frame(corr=rep(1,add.node),
                                     p.value=rep(0.05,add.node),
                                     color=add.node.color,
                                     frame.color=add.node.color)


           node.display_comb<-rbind(node.display,add.node.data)

           igraph::V(g)$color <- node.display_comb$color
           igraph::V(g)$frame.color <- node.display_comb$frame.color
         }
         ###grid节点大小
         V(g)$frame.size<-max(V(g)$size)

         newsize<-V(g)$frame.size
         newsize[diag] <-V(g)$size[diag]
         newsize[(nrow(mat2)-add.node+1):nrow(mat2)] <-0

         newcolor<-rep(grid.color,nrow(mat2))
         newcolor[diag] <-diag.frame.color
         newcolor[(nrow(mat2)-add.node+1):nrow(mat2)] <-mantel.frame.color
         #??Plot.igraph
         ###add background grid
         g1<-g

         g2<-g1%>%add_vertices(nrow(mat2),
                               frame.color=newcolor,
                               name="",
                               frame.size=newsize,
                               size=newsize,
                               weight=newsize)


         sub_net_layout2<-rbind(sub_net_layout,sub_net_layout)

         if (add.sig) {
           corp<-pp$corrPos[,1:6]
           corp$p.value<-corp$p.value%>%as.double()
           corp[which(corp$p.value<=0.001 & corp$p.value >0),]$p.value <-"***"
           corp[which(corp$p.value<=0.01 & corp$p.value >0.001),]$p.value <-"**"
           corp[which(corp$p.value<0.05 & corp$p.value >0.01),]$p.value <-"*"
           corp[which(corp$p.value>=0.05 | corp$p.value==0),]$p.value<-""
           corp<-corp$p.value%>%as.matrix()

           newp<-corp
           newp[diag] <-""
           newp[(nrow(mat2)-add.node+1):nrow(mat2)] <-""
           vertex_attr(g2)$name[seq(nrow(mat2)+1,nrow(mat2)+nrow(mat2))]<-newp


           label.dist1<-label.dist
           label.dist1[seq(1,nrow(mat2))]<-0
           label.dist1[diag] <- 0
           label.dist1[(nrow(mat2)-add.node+1):nrow(mat2)] <-0
           label.dist2<-c(label.dist,label.dist1)
         }

         layout=layout
         if (layout == "lower") {

           edge.curved1<- matrix(rep(1,add.node),ncol(netp),nrow(netp))
           for (jj in 1:add.node) {
             order <- jj*2
             edge.curved<-rep(1,nrow(corr))
             edge.curved[order]<-0
             edge.curved[1:order]<-seq(edge.curve.max,0,length.out=order)
             edge.curved[order:length(edge.curved)]<-seq(0,edge.curve.min,length.out=length(edge.curved)-order+1)
             edge.curved1[,jj]<-edge.curved
           }

           edge.curved1<-t(edge.curved1)
           E(g1)$curved<-edge.curved1

           label.locs1<-label.locs*(-1)
           label.dist1<-label.dist*(-1)
           sub_net_layout1<-sub_net_layout%>%data.frame()
           sub_net_layout1$x<-nrow(corr)+3-sub_net_layout1$x
           sub_net_layout1$y<-nrow(corr)+3-sub_net_layout1$y
           sub_net_layout1<-sub_net_layout1%>%as.matrix()

           label.dist1<-c(label.dist1,label.dist2[(nrow(mat2)+1):(nrow(mat2)+nrow(mat2))])

           ##vertex.attributes(g2)
           plot(g1,edge.curved = edge.curved1,
                vertex.shape=shapes,
                vertex.label=vertex_attr(g2)$name,
                layout=sub_net_layout1,
                vertex.label.degree = label.locs1,
                vertex.label.dist = label.dist1,
                vertex.label.family = "serif")
         }

         if (layout == "upper") {##vertex.attributes(g2)
           edge.curved1<- matrix(rep(1,add.node),ncol(netp),nrow(netp))
           for (jj in 1:add.node) {
             order <- jj*2
             edge.curved<-rep(1,nrow(corr))
             edge.curved[order]<-0
             edge.curved[1:order]<-seq(edge.curve.max,0,length.out=order)
             edge.curved[order:length(edge.curved)]<-seq(0,edge.curve.min,length.out=length(edge.curved)-order+1)
             edge.curved1[,jj]<-edge.curved
           }

           edge.curved1<-t(edge.curved1)
           E(g1)$curved<-edge.curved1

           plot(g1,edge.curved = edge.curved1,
                vertex.shape=shapes,
                vertex.label=vertex_attr(g2)$name,
                layout=sub_net_layout2,
                vertex.label.degree = label.locs,
                vertex.label.dist = label.dist2,
                vertex.label.family = "serif")}}
}


"matpMerge" <- function(corp,mantelp) {
  d1<-corp
  d5<-mantelp

  add.param<-colnames(d5)[1:length(colnames(d5))]

  add.num<-length(add.param)
  net_m<-d5%>%as.matrix()%>%abs()
  net<-d1

  add.order<-match(rownames(net),rownames(net_m))
  net_m<-net_m[add.order,]


  return(net_m)
}

"matrMerge" <- function(corp,mantelp) {
  d1<-corp
  d5<-mantelp

  add.param<-colnames(d5)[1:length(colnames(d5))]

  add.num<-length(add.param)
  net_m<-d5%>%as.matrix()%>%abs()
  net<-d1

  add.order<-match(rownames(net),rownames(net_m))
  net_m<-net_m[add.order,]


  return(net_m)
}
