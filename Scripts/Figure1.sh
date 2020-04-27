# Figure 1a
```{r}
require(ggplot2)

data <- read.table('Fig1ab.Fig2a.ED1i.ED3a.ED6d_WT_UMAP.tsv', header=T)

ggplot() + geom_point(data=subset(data, stage %in% c('WT70','WT75','WT80','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT65'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position="none")
ggplot() + geom_point(data=subset(data, stage %in% c('WT65','WT75','WT80','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT70'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position="none")
ggplot() + geom_point(data=subset(data, stage %in% c('WT65','WT70','WT80','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT75'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position="none")
ggplot() + geom_point(data=subset(data, stage %in% c('WT65','WT70','WT75','WT85')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT80'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position="none")
ggplot() + geom_point(data=subset(data, stage %in% c('WT65','WT70','WT75','WT80')), aes(x=UMAP1, y=UMAP2), color='gray85', size=0.1, alpha=0.5) + geom_point(data=subset(data, stage == 'WT85'), aes(x=UMAP1, y=UMAP2), size=0.1, alpha=0.5, color='black') + theme_void() + theme(legend.position="none")
```


# Figure 1b
```{r}
require(ggplot2)

cols <- data.frame(cluster=c(17, 4, 25, 8, 37, 3, 2, 38, 36, 11, 26, 19, 13, 20, 29, 27, 14, 0, 5, 7, 9, 28, 1, 30, 12, 23, 15, 34, 35, 33, 24, 39, 22, 31, 41, 18, 32, 16, 21, 10, 40, 6), color=c('#B9D8A5', '#FFD100', '#CBBBA0', '#92C464', '#DF76AC', '#A2C5EA', '#BD6FAB', '#477EC0', '#5B57A2', '#107F71', '#00A19A', '#3AAA35', '#7A2182', '#783F91', '#F9B233', '#E20146', '#E94E1B', '#9E937E', '#B6A9D3', '#FFED00', '#581D6F', '#6C6556', '#30892B', '#0070B2', '#954B97', '#5567AE', '#B670AC', '#8C1D82', '#2FAC66', '#1C4024', '#1A5B1A', '#277424', '#662483', '#006F72', '#6B3064', '#66296B', '#64358C', '#009767', '#642F2C', '#4A712E', '#F39200', '#624758'))

data <- read.table('Fig1ab.Fig2a.ED1i.ED3a.ED6d_WT_UMAP.tsv', header=T)
data$cluster <- factor(data$cluster, levels=cols$cluster)

ggplot(data, aes(x=UMAP1, y=UMAP2, color=cluster), size=0.1, alpha=0.5) + geom_point() + theme_void() + theme(legend.position='none') + scale_color_manual(values=as.character(cols$color))

```
