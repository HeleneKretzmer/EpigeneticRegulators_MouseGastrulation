# ED 1b
```{r}
require(RColorBrewer)

col_pals <- brewer.pal.info[brewer.pal.info$category=='qual',]
color <- unlist(mapply(brewer.pal, col_pals$maxcolors, rownames(col_pals)))

data_left <- read.table('PCA_Clustering.tsv', header=T, row.names=1)
data_center <- read.table('PCA_Unstable.tsv', header=T, row.names=1)
data_right <- read.table('PCA_Filtering.tsv', header=T, row.names=1)

plot(data_left[,1:2], col=color[data_left$embryo], main='PCA: SNP profiles of all cells')

plot(data_left[,1:2], col=color[data_left$embryo], main='PCA: SNP profiles of all cells')
points(data_center[,1:2], col=c('black','red')[data_center$embryo])

plot(data_right[,1:2], col=color[data_right$embryo], main='PCA: SNP profiles of stable cells')
```


# ED 1c
```{r}
require(ggplot2)

data <- read.table('WT_sex_expr.tsv', header=T)

ggplot(melt(data, id.vars=c('embryo', 'stage')), aes(x=embryo, y=value, fill=variable)) + geom_bar(position='dodge', stat='identity') + theme_classic() + facet_grid(~stage)
```


# ED 1e left
```{r}
data <- read.table('Marker_count_vs_unique.tsv', header=T)

plot(data$marker_count, data$fraction_unique)
```


# ED 1e right
```{r}
require(ggplot2)
require(gridExtra)

data <- read.table('Marker_count_vs_top30.tsv', header=T)

data.all <- subset(data, case=='all')
data.top30 <- subset(data, case=='top30')

data.all$cluster <- factor(data.all$cluster, levels = data.all$cluster[order(-data.all$count)])
data.top30$cluster <- factor(data.top30$cluster, levels = data.top30$cluster[order(-data.top30$count)])

p1 <- ggplot() + geom_line(data=data.all, aes(cluster, count, group=as.factor(1)), size=2) + theme_classic() + ylim(c(0,0.08))
p2 <- ggplot() + theme_classic() + geom_line(data=data.top30, aes(cluster, count, group=as.factor(1)), size=2, color='red') + ylim(c(0,0.08))
grid.arrange(p1, p2, ncol=2)
```


# ED 1f
```{r}
require(ggplot2)

data <- read.table('WTClosestCluster.tsv', header=T)

ggplot(melt(data, id.vars=c('type','cluster')), aes(x=as.factor(cluster), fill=variable, y=value)) + geom_boxplot(outlier.size=0.05) + theme_classic() + scale_fill_manual(values=c('forestgreen','grey'))
```


# ED 1g
```{r}
require(ggplot2)

cluster_colors <- c('black','#ffd100','#cbbba0','#92c464','#df76ac','#a2c5ea','#bd6fab','#477ec0','#5b57a2','#107f71','#00a19a','#3aaa35','#7a2182','#783f91','#f9b233','#e20146','#e94e1b','#9e937e','#b6a8d3','#ffed00','#581d6f','#6c6556','#30892b','#0070b2','#954b97','#5567ae','#b670ac','#8c1d82','#2fac66','#1c4024','#1a5b1a','#277424','#662483','#006f72','#6b3064','#66296b','#64358c','#009767','#642f2c','#4a712e','#f39200','#624758','#624758')
cluster_order <- c('17','4','25','8','37','3','2','38','36','11','26','19','13','20','29','27','14','0','5','7','9','28','1','30','12','23','15','34','35','33','24','39','22','31','41','18','32','16','21','10','40','6')

data <- read.table('CellStateProportions_per_Embryo.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cluster_order)

ggplot(subset(data, type %in% 'WT'), aes(x=embryo, y=frac, fill=cluster)) + geom_bar(stat='identity', size=0.5) + theme_classic(16) + facet_grid(type~stage, scales='free_y', space='free_y') + scale_fill_manual(values=cluster_colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```


# ED 1h left
```{r}
data <- read.table('CellState_prevalence.tsv', header=T)

data$cluster <- factor(data$cluster, levels=c('6','40','10','21','16','32','18','41','31','22','39','24','33','35','34','15','23','12','30','1','28','9','7','5','0','14','27','29','20','13','19','26','11','36','38','2','3','37','8','25','4','17'))

ggplot(data, aes(x=stage, y=cluster)) + geom_tile(aes(fill=norm)) + theme_classic() + scale_fill_gradient(low='white', high='black')
```


# ED 1h right
```{r}
require(Seurat)

load('WTheatmap.Robj')

DoHeatmap(data, use.scaled=T, slim.col.label=T, group.by = 'ident', group.order = rev(c.order), draw.line = FALSE, cex.row=1, group.spacing=0.05, group.cex=5, col.low='dodgerblue1', col.mid='black', col.high='yellow', genes.use=genes)

```


# ED 1i left
```{r}
require(ggplot2)

data <- read.table('WT_UMAP.tsv', header=T)

ggplot(data, aes(x=UMAP1, y=UMAP2, color=stage)) + geom_point(size=0.1) + theme_void() + theme(legend.position='none') + scale_color_grey()
```


# ED 1i right
```{python}
import scanpy as sc
import scvelo as scv
scv.settings.set_figure_params('scvelo')

WT = sc.read('WT_velocity.h5ad', sparse=True)
adata.obs['CellType'].cat.set_categories(['17','4','25','8','37','3','2','38','36','11','26','19','13','20','29','27','14','0','5','7','9','28','1','30','12','23','15','34','35','33','24','39','22','31','41','18','32','16','21','10','40','6'], inplace=True)

scv.pl.velocity_embedding_stream(WT, legend_loc='on data', alpha=0.5, color='CellType')
```
