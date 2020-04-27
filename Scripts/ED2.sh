# ED 2b
```{r}
require(ggplot2)

data <- read.table('ED2b_KO_expression_in_WT_per_stage.tsv', header=T)

ggplot(data, aes(x=stage, y=frac, fill=stage)) + geom_boxplot(outlier.shape=NA) + geom_point(size=0.5) + theme_classic() + facet_grid(~type.y) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_grey(start = 0.8, end = 0.1)
```


# ED 2c
```{r}
require(ggplot2)

data <- read.table('ED2c_KO_expression_in_WT_per_lineage.tsv', header=T)

ggplot(data, aes(x=lineage, y=value, fill=lineage)) + geom_boxplot() + theme_classic() + facet_grid(~variable) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c('black','#E6007E', '#5BAE6F','#FFED45','#5F2C81', '#C98DA5','#CA9E67','#B090BF'))
```


# ED 2d
```{r}
require(ggplot2)

data <- read.table('ED2d_DNMT1_85_statistics.tsv', header=T)
data <- read.table('ED2d_DNMT3A_85_statistics.tsv', header=T)
data <- read.table('ED2d_DNMT3B_85_statistics.tsv', header=T)
data <- read.table('ED2d_G9A_85_statistics.tsv', header=T)
data <- read.table('ED2d_EED_85_statistics.tsv', header=T)
data <- read.table('ED2d_RNF2_85_statistics.tsv', header=T)
data <- read.table('ED2d_KDM2B_85_statistics.tsv', header=T)
data <- read.table('ED2d_L3MBTL2_85_statistics.tsv', header=T)
data <- read.table('ED2d_KMT2A_85_statistics.tsv', header=T)
data <- read.table('ED2d_KMT2B_85_statistics.tsv', header=T)

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


# ED 2e
```{r}
require(ggplot2)

data <- read.table('ED2e_ClosestCluster.tsv', header=T)

ggplot(melt(data), aes(x=type, fill=variable, y=value)) + geom_boxplot(outlier.size=0.05) + theme_classic() + scale_fill_manual(values=c('forestgreen','grey'))
```


# ED 2f
```{r}
require(ggplot2)

cluster_colors <- c('black','#ffd100','#cbbba0','#92c464','#df76ac','#a2c5ea','#bd6fab','#477ec0','#5b57a2','#107f71','#00a19a','#3aaa35','#7a2182','#783f91','#f9b233','#e20146','#e94e1b','#9e937e','#b6a8d3','#ffed00','#581d6f','#6c6556','#30892b','#0070b2','#954b97','#5567ae','#b670ac','#8c1d82','#2fac66','#1c4024','#1a5b1a','#277424','#662483','#006f72','#6b3064','#66296b','#64358c','#009767','#642f2c','#4a712e','#f39200','#624758','#624758')
cluster_order <- c('17','4','25','8','37','3','2','38','36','11','26','19','13','20','29','27','14','0','5','7','9','28','1','30','12','23','15','34','35','33','24','39','22','31','41','18','32','16','21','10','40','6')

data <- read.table('ED1g.ED2f_CellStateProportions_per_Embryo.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cluster_order)

ggplot(subset(data, stage %in% '85'), aes(x=embryo, y=frac, fill=cluster)) + geom_bar(stat='identity', size=0.5) + theme_classic(16) + facet_grid(type~stage, scales='free_y', space='free_y') + scale_fill_manual(values=cluster_colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```
