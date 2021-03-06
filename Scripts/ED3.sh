# ED 3a
```{r}
require(ggplot2)

data <- read.table('Fig1ab.Fig2a.ED1i.ED3a.ED6d_WT_UMAP.tsv', header=T)

ggplot() + geom_point(data=subset(data, stage %in% c('WT70','WT75','WT80','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT65'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position='none')
ggplot() + geom_point(data=subset(data, stage %in% c('WT65','WT75','WT80','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT70'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position='none')
ggplot() + geom_point(data=subset(data, stage %in% c('WT65','WT70','WT80','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT75'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position='none')
ggplot() + geom_point(data=subset(data, stage %in% c('WT65','WT70','WT75','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT80'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position='none')
ggplot() + geom_point(data=subset(data, stage %in% c('WT65','WT70','WT75','WT80')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT85'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position='none')


cluster_colors <- c('black','#ffd100','#cbbba0','#92c464','#df76ac','#a2c5ea','#bd6fab','#477ec0','#5b57a2','#107f71','#00a19a','#3aaa35','#7a2182','#783f91','#f9b233','#e20146','#e94e1b','#9e937e','#b6a8d3','#ffed00','#581d6f','#6c6556','#30892b','#0070b2','#954b97','#5567ae','#b670ac','#8c1d82','#2fac66','#1c4024','#1a5b1a','#277424','#662483','#006f72','#6b3064','#66296b','#64358c','#009767','#642f2c','#4a712e','#f39200','#624758','#624758')
cluster_order <- c('17','4','25','8','37','3','2','38','36','11','26','19','13','20','29','27','14','0','5','7','9','28','1','30','12','23','15','34','35','33','24','39','22','31','41','18','32','16','21','10','40','6')
s <- c(17, 27, 8,19,1,11,16,10,24,39,33,35,26,31, 3,7,38,23, 37,2,13,20,12,30,36,32,34,21,6,18,9,22, 5,15,41, 0,25,28, 4,29,40,14)

data <- read.table('Fig2f.ED3a.ED6d_MedianCellStateProportions.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cluster_order[match(s, cluster_order)])

ggplot(subset(data, !type %in% c('Cdkn2a','Eed.Cdkn2a','WT_65','WT_70','WT_75','WT_80')), aes(x=type, y=frac, fill=cluster)) + geom_bar(stat='identity', size=0.5) + theme_classic(16) + scale_fill_manual(values=cluster_colors[match(s, cluster_order)])


cols <- data.frame(cluster=c(17, 4, 25, 8, 37, 3, 2, 38, 36, 11, 26, 19, 13, 20, 29, 27, 14, 0, 5, 7, 9, 28, 1, 30, 12, 23, 15, 34, 35, 33, 24, 39, 22, 31, 41, 18, 32, 16, 21, 10, 40, 6), color=c('#B9D8A5', '#FFD100', '#CBBBA0', '#92C464', '#DF76AC', '#A2C5EA', '#BD6FAB', '#477EC0', '#5B57A2', '#107F71', '#00A19A', '#3AAA35', '#7A2182', '#783F91', '#F9B233', '#E20146', '#E94E1B', '#9E937E', '#B6A9D3', '#FFED00', '#581D6F', '#6C6556', '#30892B', '#0070B2', '#954B97', '#5567AE', '#B670AC', '#8C1D82', '#2FAC66', '#1C4024', '#1A5B1A', '#277424', '#662483', '#006F72', '#6B3064', '#66296B', '#64358C', '#009767', '#642F2C', '#4A712E', '#F39200', '#624758'))

data <- read.table('Fig2a.ED3a.ED6d_KO_UMAP.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cols$cluster)

ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Dnmt3a'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Dnmt3a')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Dnmt3b'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Dnmt3b')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Dnmt1'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Dnmt1')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Eed'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Eed')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Rnf2'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Rnf2')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Kdm2b'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Kdm2b')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'L3mbtl2'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'L3mbtl2')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'G9a'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'G9a')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Kmt2a'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Kmt2a')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Kmt2b'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Kmt2b')$cluster)$color))
```

# ED 3b
```{r}
require(ggplot2)

data <- read.table('ED3b_PCA_embryo_WT.tsv', header=T)

ggplot(data, aes(x=stage, y=PC1.binary)) + geom_boxplot(outlier.shape=NA) + geom_point() + theme_classic()
ggplot(data, aes(x=stage, y=PC1.prop)) + geom_boxplot(outlier.shape=NA) + geom_point() + theme_classic()

```


# ED 3c
```{r}
require(ggplot2)

data <- read.table('ED3c.ED3d_PCA_median_WT.tsv', header=T)

ggplot(pca[c(1:5),], aes(x=PC1.binary, y=PC1.prop)) + geom_point() + theme_classic()

```


# ED 3d
```{r}
require(ggplot2)

dataWT <- read.table('ED3c.ED3d_PCA_median_WT.tsv', header=T)
dataKO <- read.table('ED3d_PCA_embryo_KO.tsv', header=T)

ggplot() + geom_point(data=dataWT[,-1], aes(x=PC1.binary, y=PC1.prop), shape=8) + theme_classic() + geom_point(data=dataKO, aes(x=PC1.binary, y=PC1.prop, group=type), size=2, shape=22, fill='firebrick') + facet_grid(~type)

```


# ED 3e top
```{r}
require(ComplexHeatmap)
require(circlize)
require(cluster)

data <- read.table('Fig2c.ED3e_CellStateFraction_diff_expr_E.tsv', header=T)

col_fun_2 = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c('blue','blue', 'grey95', 'red','red'))

Heatmap(as.matrix(data), show_row_names=FALSE, name='Fraction embryonic\ncell states', col=col_fun_2, cluster_columns=as.dendrogram(agnes(t(as.matrix(data)))))
```


# ED 3e bottom
```{r}
require(ComplexHeatmap)
require(circlize)
require(cluster)

data <- read.table('ED3e_CellStateFraction_diff_expr_X.tsv', header=T)

col_fun_2 = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c('blue','blue', 'grey95', 'red','red'))

Heatmap(as.matrix(data), show_row_names=FALSE, name='Fraction extraembryonic\ncell states', col=col_fun_2, cluster_columns=as.dendrogram(agnes(t(as.matrix(data)))))
```






