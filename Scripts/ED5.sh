# ED 5a top
```{r}
require(ggplot2)
require(RColorBrewer)

data <- read.table('RepeatFamily_methylation.tsv', header=T)

ggplot(melt(data), aes(x=variable, y=value)) + geom_point(aes(color=Family), shape=4, size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90)) + facet_grid(~gsub('.*_','',variable), scales='free_x', space='free_x') + ylim(c(0,1)) + scale_color_manual(values=c(brewer.pal(6,'Blues')[-c(1,2)], brewer.pal(3,'Reds')[-1]))
```


# ED 5a bottom
```{r}
require(ggplot2)
require(RColorBrewer)

data <- read.table('RepeatFamily_expression.tsv', header=T)

ggplot(data, aes(x=type, y=avg_expr, color=element)) + geom_point(shape=4, size=2) + theme_classic() + facet_grid(~tissue) + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylim(c(0,0.015)) + scale_color_manual(values=c(brewer.pal(6,'Blues')[-c(1,2)], brewer.pal(3,'Reds')[-1]))
```


# ED 5b top
```{r}
require(vioplot)

data <- read.table('IAPEz-int_methylation.tsv', header=T)

vioplot(data[, grep('Epi', colnames(data))], names=c(colnames(data[, grep('Epi', colnames(data))],)), cex=0.5, main='IAPEz-int')
vioplot(data[, grep('Exe', colnames(data))], names=c(colnames(data[, grep('Exe', colnames(data))],)), cex=0.5, main='IAPEz-int')
```


# ED 5b bottom
```{r}
require(ggplot2)

data <- read.table('IAPEz-int_expression.tsv', header=T)

ggplot(data, aes(x=type, y=fraction)) + geom_point(shape=4, size=2) + theme_classic() + facet_grid(~tissue) + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylim(c(0,0.015))
```


# ED 5c
```{r}
require(ggplot2)

data <- read.table('L3mbtl2_methyl_expr.tsv', header=T)

ggplot() +
geom_hline(yintercept=c(-0.2,0.2), color='grey') + geom_vline(xintercept=c(-0.1,0.1), color='grey') + theme_classic() +
geom_point(data=subset(data, XE=="X" & abs(delta_Epi)<0.1 & abs(diff)<0.2), aes(x=delta_Epi, y=diff), color='grey', size=1.5) +
geom_point(data=subset(data, XE=="X" & (abs(delta_Epi)>=0.1 | abs(diff)>=0.2)), aes(x=delta_Epi, y=diff), color='black', size=2) +
geom_point(data=subset(data, XE=="X" & delta_Epi <= -0.1 & diff >= 0.2), aes(x=delta_Epi, y=diff), color='forestgreen', size=3) +
geom_text(data=subset(data, XE=="X" & delta_Epi <= -0.1 & diff >= 0.2), aes(x=delta_Epi+0.05, y=diff, label=gene), color='forestgreen', size=2) + ggtitle('Extraembryonic') + xlim(c(-1,1)) + ylim(c(-1,1)) + xlab('delta_methyl') + ylab('delta_expr')
```


# ED 5e top
```{r}
require(vioplot)

data <- read.table('L3mbtl2_hypo_derepr_genes_methylation.tsv', header=T, row.names=1)

vioplot(data, names=colnames(data))
```


# ED 5e bottom
```{r}
require(ggplot2)

data <- read.table('L3mbtl2_hypo_derepr_genes_pos_cells.tsv', header=T)

ggplot(data, aes(x=KO, y=pos.cells)) + geom_boxplot() + theme_classic() + facet_wrap(~XE) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```
