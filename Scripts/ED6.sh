# ED 6a
```{r}
require(ggplot2)
require(reshape2)

data <- read.table('ED6a_PCG_dist_to_CellState.tsv', header=T)

ggplot(melt(data), aes(x=type, fill=variable, y=value)) + geom_boxplot(outlier.size=0.05) + theme_classic() + scale_fill_manual(values=c('forestgreen','grey','darkgrey'))
```


# ED 6b
```{r}
require(ggplot2)
require(reshape2)

data <- read.table('Fig3e.ED6b.ED10f_X_expr.tsv', header=T)

ggplot(subset(data, stage %in% c('WT_65','WT_70','WT_75','WT_80','WT_85','Eed_85','Rnf2_85','Kdm2b_85')), aes(x=gsub('_.*','',stage), y=log2(frac), fill=sex)) + geom_boxplot(width=0.5, outlier.shape=NA) + theme_classic() + scale_fill_manual(values=c('firebrick','dodgerblue')) + facet_grid(~XE)
```


# ED 6c
```{r}
require(ggplot2)

data <- read.table('ED6c_CDKN2A_85_Cdkn2a_statistics.tsv', header=T)
data <- read.table('ED6c_EED.CDKN2A_85_Cdkn2a_statistics.tsv', header=T)
data <- read.table('ED6c_EED.CDKN2A_85_Eed_statistics.tsv', header=T)

if('g4' %in% data$g){
	data$ID <- factor(data$ID, levels=c(EedWT_g1Eed, EedKO_g1Eed, EedWT_g2Eed, EedKO_g2Eed, EedWT_g3Eed, EedKO_g3Eed, 'WT_g4', 'KO_g4'))
	data$read <- factor(data$read, levels=c('Mismatched', 'Spliced', 'Insufficient', 'Complete'))
	ggplot(df, aes(x=ID, y=Freq, fill=read)) + geom_bar(stat='identity', width=0.5, color='black') + theme_classic() + facet_grid(~g, space='free_x', scales='free_x') + scale_fill_manual(values=rev(c('lightblue','lightgrey','darkgrey','darkblue'))) + xlab('') + ylab('Number of reads')
}else{
	data$ID <- factor(data$ID, levels=c('WT_g1', 'KO_g1', 'WT_g2', 'KO_g2', 'WT_g3', 'KO_g3'))
	data$read <- factor(data$read, levels=c('Mismatched', 'Spliced', 'Insufficient', 'Complete'))
	ggplot(data, aes(x=ID, y=Freq, fill=read)) + geom_bar(stat='identity', width=0.5, color='black') + theme_classic() + facet_grid(~g, space='free_x', scales='free_x') + scale_fill_manual(values=rev(c('lightblue','lightgrey','darkgrey','darkblue'))) + xlab('') + ylab('Number of reads')
}
```




# ED 6d
```{r}
require(ggplot2)

cols <- data.frame(cluster=c(17, 4, 25, 8, 37, 3, 2, 38, 36, 11, 26, 19, 13, 20, 29, 27, 14, 0, 5, 7, 9, 28, 1, 30, 12, 23, 15, 34, 35, 33, 24, 39, 22, 31, 41, 18, 32, 16, 21, 10, 40, 6), color=c('#B9D8A5', '#FFD100', '#CBBBA0', '#92C464', '#DF76AC', '#A2C5EA', '#BD6FAB', '#477EC0', '#5B57A2', '#107F71', '#00A19A', '#3AAA35', '#7A2182', '#783F91', '#F9B233', '#E20146', '#E94E1B', '#9E937E', '#B6A9D3', '#FFED00', '#581D6F', '#6C6556', '#30892B', '#0070B2', '#954B97', '#5567AE', '#B670AC', '#8C1D82', '#2FAC66', '#1C4024', '#1A5B1A', '#277424', '#662483', '#006F72', '#6B3064', '#66296B', '#64358C', '#009767', '#642F2C', '#4A712E', '#F39200', '#624758'))


data <- read.table('Fig1ab.Fig2a.ED1i.ED3a.ED6d_WT_UMAP.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cols$cluster)

ggplot() + geom_point(data=subset(data, stage %in% c('WT65', 'WT70','WT75','WT80','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT85'), aes(x=UMAP1, y=UMAP2, color=as.factor(cluster)), size=0.1, alpha=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% unique(subset(data, stage == 'WT85'))$cluster)$color))


data <- read.table('Fig2a.ED3a.ED6d_KO_UMAP.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cols$cluster)

ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Cdkn2a'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Cdkn2a')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Eed.Cdkn2a'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Eed.Cdkn2a')$cluster)$color))
ggplot() + geom_point(data=subset(data, type =='WT'), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, type == 'Eed'), aes(x=UMAP1, y=UMAP2, color=cluster), size=0.5) + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(subset(cols, cluster %in% subset(data, type == 'Eed')$cluster)$color))


cluster_colors <- c('black','#ffd100','#cbbba0','#92c464','#df76ac','#a2c5ea','#bd6fab','#477ec0','#5b57a2','#107f71','#00a19a','#3aaa35','#7a2182','#783f91','#f9b233','#e20146','#e94e1b','#9e937e','#b6a8d3','#ffed00','#581d6f','#6c6556','#30892b','#0070b2','#954b97','#5567ae','#b670ac','#8c1d82','#2fac66','#1c4024','#1a5b1a','#277424','#662483','#006f72','#6b3064','#66296b','#64358c','#009767','#642f2c','#4a712e','#f39200','#624758','#624758')
cluster_order <- c('17','4','25','8','37','3','2','38','36','11','26','19','13','20','29','27','14','0','5','7','9','28','1','30','12','23','15','34','35','33','24','39','22','31','41','18','32','16','21','10','40','6')
s <- c(17, 27, 8,19,1,11,16,10,24,39,33,35,26,31, 3,7,38,23, 37,2,13,20,12,30,36,32,34,21,6,18,9,22, 5,15,41, 0,25,28, 4,29,40,14)

data <- read.table('MedianCellStateProportions.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cluster_order[match(s, cluster_order)])

ggplot(subset(data, type %in% c('WT_85','Eed','Cdkn2a','Eed.Cdkn2a')), aes(x=type, y=frac, fill=cluster)) + geom_bar(stat='identity', size=0.5) + theme_classic(16) + scale_fill_manual(values=cluster_colors[match(s, cluster_order)])
```


# ED 6e
```{r}
require(pheatmap)

data <- read.table('ED6e_CellStateProportion_correlation_Eed_Cdkn2a.tsv', header=T, row.names=1)

pheatmap(data, cutree_rows = 3)
```

