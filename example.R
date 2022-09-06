library(mantelVisual)
library(igraph)
data(edmat)
data(param)
data("data_node")

data_test<-cor_mantel(edmat=edmat,ctmat=param,filter=FALSE,ed.p=0.05,ed.r=0.3,
                      method="pearson",
                      permutations=999,
                      mantel_select=list(afirm = 2:2,cya= 1:1,
                                         actino=4:4,
                                         bacter=5:5,
                                         pro=3:3))
data12<-data_test

{'Error in cor(as.vector(xdis), ydis, method = method, use = use) :
  cov/cor中有遗漏值'} ##数据行和为0

plot_mantel_heatmap(data12,
                    add.sig=TRUE,  ###是否添加显著性标志
                    edge.curve.min = -0.2,
                    edge.curve.max = -0.2, ###线条曲率正负阈值范围
                    layout = "lower",     ###上三角还是下三角
                    grid.background = FALSE, ###添加最值网格
                    grid.color = "grey90",  ###grid颜色
                    space.h=1.6,
                    space.v=1.5,

                    adjust.diag=TRUE, ###是否调整对角线距离
                    adjust.diag.para=TRUE, ###是否平行调整对角线距离
                    adjust.diag.para.dist=1,  ###是否平行调整对角线距离为1
                    adjust.node = c(9),   ###adjust.diag.para=false后，调整对角线节点坐标，默认为横坐标，
                    ###可以是数值也可是向量
                    adjust.diag.dist = c(0.2),  ###调整的距离
                    diag.node.size = 5,         ###对角线节点大小
                    diag.node.shape.order = 9,  ###对角线节点形状
                    diag.frame.color = "red",  #对角线边框颜色
                    diag.node.color = "red",

                    mantel.node.pos.para=TRUE, ###false,不平行，则需要调整下面参数
                    mantel.node.x.pos=c(5,7,9,11),
                    mantel.node.y.pos=c(13,11,9,7),  ###调整节点坐标，可以是数值也可是向量
                    mantel.node.size.para=TRUE,     ###mantel节点大小是否调整一致
                    mantel.node.size=5,     ###mantel.node.size.para=TRUE,mantel节点大小是否调整一致
                    mantel.node.size.custom=c(10,1,2,5,8),  ###mantel.node.size.para=FALSE5为大小，4为数量
                    mantel.node.shape.order = c(1,1,1,1,1), ###数值or向量，mantel节点形状，看shapes()函数
                    mantel.frame.color="purple", #mantel节点边框颜色
                    mantel.node.color = "purple",

                    heatmap.node.size.ratio=15,
                    heatmap.shape = 3,    ###热图形状
                    node.color.postive = "purple", ###热图正负相关颜色
                    node.color.mid = "white",
                    node.color.negative = "darkgreen",

                    adjust.text.x = 0.7,  ###控制文字距离
                    adjust.text.y = 0.7,

                    #label.cex=3,
                    #edge.size.max = 6,  ###定义边的粗细范围
                    #edge.size.center= 3,
                    #edge.size.min = 1,

                    pbreaks = c(-Inf,0.001, 0.01, 0.05, Inf),
                    plabels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05"),
                    color.p=c("darkgreen","purple","palegreen","grey90"),
                    #adjustcolor("white",alpha.f = 0.5)
                    rbreaks = c(-Inf,0, 0.2, 0.4, Inf),
                    rlabels = c("< 0","0 - 0.2", "0.2 - 0.4", ">= 0.4"),
                    size.r=c(0.5,1,2,4),
                    node.label.color=NULL)

plot_mantel_1hierarchy(data_test,
                       data_node, ###计算节点大小,0-1标准化的数据
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
                       color.p=c("darkgreen","purple","palegreen","grey90"),
                       #adjustcolor("white",alpha.f = 0.5)
                       rbreaks = c(-Inf,0, 0.2, 0.4, Inf),
                       rlabels = c("< 0","0 - 0.2", "0.2 - 0.4", ">= 0.4"),
                       size.r=c(0.5,1,2,4),
                       node.label.color=NULL,
                       vertex.receiver = seq(1,6),
                       color.use=colorCustom(50,pal = "gygn"))

par(mfrow=c(1,2))
corrplot::corrplot(corr=data_test$occor.r,p.mat=data_test$occor.p,
                   mar = c(3, 0, 0, 0),
                   method = "square",type = "full")

mantelVisual::plot_mantel_2hierarchy(data_test,
                       data_node, ###计算节点大小
                       margin=-0.2,
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
                       color.p=c("darkgreen","purple","palegreen","grey90"),
                       #adjustcolor("white",alpha.f = 0.5)
                       rbreaks = c(-Inf,0, 0.2, 0.4, Inf),
                       rlabels = c("< 0","0 - 0.2", "0.2 - 0.4", ">= 0.4"),
                       size.r=c(0.5,1,2,4),
                       node.label.color=NULL,
                       vertex.receiver = seq(1,6),
                       color.use=colorCustom(50,pal = "gygn"))
