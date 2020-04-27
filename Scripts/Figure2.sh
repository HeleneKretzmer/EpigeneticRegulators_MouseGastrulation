# Figure 2a
```{r}
require(ggplot2)

cols <- data.frame(cluster=c(17, 4, 25, 8, 37, 3, 2, 38, 36, 11, 26, 19, 13, 20, 29, 27, 14, 0, 5, 7, 9, 28, 1, 30, 12, 23, 15, 34, 35, 33, 24, 39, 22, 31, 41, 18, 32, 16, 21, 10, 40, 6), color=c('#B9D8A5', '#FFD100', '#CBBBA0', '#92C464', '#DF76AC', '#A2C5EA', '#BD6FAB', '#477EC0', '#5B57A2', '#107F71', '#00A19A', '#3AAA35', '#7A2182', '#783F91', '#F9B233', '#E20146', '#E94E1B', '#9E937E', '#B6A9D3', '#FFED00', '#581D6F', '#6C6556', '#30892B', '#0070B2', '#954B97', '#5567AE', '#B670AC', '#8C1D82', '#2FAC66', '#1C4024', '#1A5B1A', '#277424', '#662483', '#006F72', '#6B3064', '#66296B', '#64358C', '#009767', '#642F2C', '#4A712E', '#F39200', '#624758'))


data <- read.table('Fig1ab.Fig2a.ED1i.ED3a.ED6d_WT_UMAP.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cols$cluster)

ggplot() + geom_point(data=subset(data, stage %in% c('WT65', 'WT70','WT75','WT80','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT85'), aes(x=UMAP1, y=UMAP2, color=as.factor(cluster)), size=0.1, alpha=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% unique(subset(data, stage == 'WT85'))$cluster)$color))


data <- read.table('KO_UMAP.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cols$cluster)

ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Eed'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Eed')$cluster)$color))
```


# Figure 2b
```{r}
require(scatterpie)
require(ggplot2)

data <- read.table('Fig2b_KO_embryo_staging.tsv', header=T)

ggplot() + geom_scatterpie(data=data, aes(x=x, y=y, r=log(1.1+r)), cols=c('female','male')) + coord_equal() + scale_x_continuous(name='KO', breaks=data$x, labels=data$type) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```


# Figure 2c top
```{r}
require(cluster)
require(circlize)
require(ComplexHeatmap)

data <- read.table('Fig2c_CellState_FC.tsv', header=T, row.names=1)

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c('blue', 'white', 'red'))
ha <- HeatmapAnnotation(foo = colMeans(data, na.rm=T), col = list(foo = col_fun))

Heatmap(data, cluster_columns=as.dendrogram(agnes(t(as.matrix(data)))), top_annotation = ha)
```


# Figure 2c bottom
```{r}
require(ComplexHeatmap)
require(circlize)
require(cluster)

data <- read.table('Fig2c_CellStateFraction_diff_expr_E.tsv', header=T)

col_fun_2 = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c('blue','blue', 'grey95', 'red','red'))

Heatmap(as.matrix(data), show_row_names=FALSE, name='Fraction embryonic\ncell states', col=col_fun_2, cluster_columns=as.dendrogram(agnes(t(as.matrix(data)))))
```



# Figure 2d
```{r}
require(ggplot2)
require(gridExtra)

CGI_Epi <- read.table('Fig2d_CGI_Epi.tsv', header=T)
CGI_Exe <- read.table('Fig2d_CGI_ExE.tsv', header=T)

p1 <- ggplot(na.omit(melt(CGI_Epi[,-c(1:3)], id.vars='genes')), aes(x=variable, fill=value)) + geom_bar() + theme_classic() + facet_grid(~gsub('_.*','',value)) + ylim(c(0,2500)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
p2 <- ggplot(na.omit(melt(CGI_Exe[,-c(1:3)], id.vars='genes')), aes(x=variable, fill=value)) + geom_bar() + theme_classic() + facet_grid(~gsub('_.*','',value)) + ylim(c(0,2500)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
grid.arrange(p1,p2,nrow=1) 
```


# Figure 2e methylation
```{r}
require(vioplot)

data <- read.table('Fig2e.ED5b_IAPEz-int_methylation.tsv', header=T)

vioplot(data[, grep('Epi', colnames(data))], names=c(colnames(data[, grep('Epi', colnames(data))],)), cex=0.5, main='IAPEz-int')
vioplot(data[, grep('Exe', colnames(data))], names=c(colnames(data[, grep('Exe', colnames(data))],)), cex=0.5, main='IAPEz-int')
```


# Figure 2e expression
```{r}
require(vioplot)

data <- read.table('Fig2e.ED5b_IAPEz-int_expression.tsv', header=T)

ggplot(subset(data, type %in% c('WT','Dnmt1','Dnmt3a','Dnmt3b','G9a')), aes(x=type, y=fraction)) + geom_point(shape=4, size=2) + theme_classic() + facet_grid(~tissue) + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylim(c(0,0.015))
```


# Figure 2f
```{r}
require(ggplot2)

cluster_colors <- c('black','burlywood1','blanchedalmond','black','black','black','black','black','black','black','black','black','black','black','burlywood1','black','burlywood1','blanchedalmond','black','black','black','blanchedalmond','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','burlywood1','black')
cluster_order <- c('17','4','25','8','37','3','2','38','36','11','26','19','13','20','29','27','14','0','5','7','9','28','1','30','12','23','15','34','35','33','24','39','22','31','41','18','32','16','21','10','40','6')
s <- c(17, 27, 8,19,1,11,16,10,24,39,33,35,26,31, 3,7,38,23, 37,2,13,20,12,30,36,32,34,21,6,18,9,22, 5,15,41, 0,25,28, 4,29,40,14)

data <- read.table('Fig2f.ED3a.ED6dMedianCellStateProportions.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cluster_order[match(s, cluster_order)])

ggplot(subset(data, type %in% c('L3mbtl2')), aes(x=type, y=frac, fill=cluster)) + geom_bar(stat='identity', size=0.5) + theme_classic(16) + scale_fill_manual(values=cluster_colors[match(s, cluster_order)]) + coord_polar('y')
ggplot(subset(data, type %in% c('WT_70')), aes(x=type, y=frac, fill=cluster)) + geom_bar(stat='identity', size=0.5) + theme_classic(16) + scale_fill_manual(values=cluster_colors[match(s, cluster_order)]) + coord_polar('y')
```


# Figure 2g
```{r}
require(ggplot2)

data <- read.table('Fig2g.ED5c_L3mbtl2_methyl_expr.tsv', header=T)

ggplot() +
geom_hline(yintercept=c(-0.2,0.2), color='grey') + geom_vline(xintercept=c(-0.1,0.1), color='grey') + theme_classic() +
geom_point(data=subset(data, XE=='E' & abs(delta_Epi)<=0.1 & abs(diff)<=0.2), aes(x=delta_Epi, y=diff), color='grey', size=1.5) +
geom_point(data=subset(data, XE=='E' & (abs(delta_Epi)>=0.1 | abs(diff)>=0.2)), aes(x=delta_Epi, y=diff), color='black', size=2) +
geom_point(data=subset(data, XE=='E' & delta_Epi <= -0.1 & diff >= 0.2), aes(x=delta_Epi, y=diff), color='forestgreen', size=3) +
geom_text(data=subset(data, XE=='E' & delta_Epi <= -0.1 & diff >= 0.2), aes(x=delta_Epi+0.05, y=diff, label=gene), color='forestgreen', size=2) + ggtitle('Embryonic') + xlim(c(-1,1)) + ylim(c(-1,1)) + xlab('delta_methyl') + ylab('delta_expr')
```

