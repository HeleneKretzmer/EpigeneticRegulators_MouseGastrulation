# ED 7a left
```{r}
require(cluster)
require(circlize)
require(ComplexHeatmap)

data_me <- read.table('DMV_methylation.txv',  header=T, row.names=1)
data_enr <- read.table('DMV_enrichment.txv',  header=T, row.names=1)


Heatmap(as.matrix(data_me), col = colorRamp2(c(0, 0.15, 0.5, 1), c('white', 'white', 'red','red')), rect_gp=gpar(col='grey', lwd=1))
Heatmap(as.matrix(data_enr), col = colorRamp2(c(0, 1, 9), c('lightgrey', 'lightgrey', 'forestgreen')), rect_gp=gpar(col='grey', lwd=1))
```


# ED 7a right
```{r}
require(vioplot)

DMV <- read.table('DMV_methyl.tsv', header=T)
CGI <- read.table('CGI_methyl.tsv', header=T)
CGIhyper <- read.table('CGIhyper_methyl.tsv', header=T)

par(mfrow=c(2,3))
vioplot(DMV[,1:4], names=colnames(DMV[,1:4]), ylim=c(0,1))
vioplot(CGI[,1:4], names=colnames(CGI[,1:4]), ylim=c(0,1))
vioplot(CGIhyper[,1:4], names=colnames(CGIhyper[,1:4]), ylim=c(0,1))
vioplot(DMV[,5:8], names=colnames(DMV[,5:8]), ylim=c(0,1))
vioplot(CGI[,5:8], names=colnames(CGI[,5:8]), ylim=c(0,1))
vioplot(CGIhyper[,5:8], names=colnames(CGIhyper[,5:8]), ylim=c(0,1))
```



# ED 7b
```{r}
require(ComplexHeatmap)
require(cluster)
require(circlize)

annotation <- data.frame(
cluster=c(17, 27, 4, 14, 29, 40, 25, 0, 28, 5, 15, 41, 37, 2, 13, 32, 36, 9, 22, 30, 26, 31, 20, 21, 6, 34, 18, 12, 3, 38, 23, 7, 8, 19, 16, 10, 11, 39, 1, 24, 33, 35),
label=c('Epi', 'PGCs', 'Xendo', 'Xendo', 'Xendo', 'Xendo', 'Xecto', 'Xecto', 'Xecto', 'Xmeso', 'Xmeso', 'Xmeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Emeso', 'Eendo', 'Eecto', 'Eecto', 'Eecto', 'Eecto', 'Eecto', 'Eecto', 'Eecto', 'Eecto', 'Eecto', 'Eecto'))

methyl <- read.table('DMV_methylation.tsv', header=T)
ht_methyl <- Heatmap(as.matrix(methyl[,c(4,6,5,7)]), cluster_rows = diana(as.matrix(methyl[,c(4,6,5,7)])), cluster_columns=F, col=colorRamp2(c(0,0.15,0.5,1), c('grey80', 'white', 'red', 'red')), show_row_dend=F, show_row_names=F)
genes <- methyl$gene[row_order(ht_methyl)]
g <- na.omit(unique(genes))

delta_expr <- read.table('DMV_PRC_recur_expr.tsv', header=T)
ht_diff <- Heatmap(delta_expr[,-1], cluster_rows=F, cluster_columns=F, show_row_names=F, col=colorRamp2(c(0,27), c('grey95', 'black')))

WT.g <- read.table('DMV_WT_expr.tsv', header=T, row.names=1)
WT.g <- WT.g[match(rownames(WT.g), genes),]
colnames(WT.g) <- gsub('X','',colnames(WT.g))

timing <- data.frame(time=c(rep('early',length(c(17,4,25,8,37,3))), rep('mid',length(c(2,38,36,11,26,19,13,20,29,27,14,0,5,7,9,28,1,30,12,23,15))),rep('late',length(c(34,35,33,24,39,22,31,41,18,32,16,21,10,40,6)))), cluster=c(c(17,4,25,8,37,3),c(2,38,36,11,26,19,13,20,29,27,14,0,5,7,9,28,1,30,12,23,15),c(34,35,33,24,39,22,31,41,18,32,16,21,10,40,6)))

ha <- HeatmapAnnotation(
    Lineage = annotation$label[match(colnames(WT.g), annotation$cluster)],
    Cluster = as.numeric(colnames(WT.g)),
    Timing = timing[match(colnames(WT.g), timing$cluster), ]$time,
    Eed = as.numeric(colnames(WT.g) %in% c(0,1,10,11,12,13,14,15,16,17,19,2,20,23,25,26,27,28,29,3,32,36,37,38,4,40,5,6,7,8,9)),
    Rnf2 = as.numeric(colnames(WT.g) %in% c(0,1,10,11,12,13,14,15,16,18,19,2,20,23,24,25,26,27,28,29,3,30,31,32,33,35,36,37,39,4,40,41,5,6,7,8,9)),
    Kdm2b = as.numeric(colnames(WT.g) %in% c(28,26,9,30,15,20,1,2,5,40,10,31,21,11,41,29,0,19,25,39,16,18,4,23,12,8,36,6,14,7,13,24,33,34,32,35,37,22,3,27,38)),
    expr_genes_freq = anno_barplot(c(26,23,36,21,55,45,13,44,51,5,38,10,37,32,26,7,27,32,21,34,30,15,19,17,45,38,34,21,30,40,9,24,19,43,22,19,48,58,45,42,3,5)),
    col = list(Lineage = c('Epi' = 'black', 'PGCs' = '#E6007E', Eecto='#5BAE6F', Eendo='#FFED45', Emeso='#5F2C81', Xecto='#C9BDA5', Xendo='#CA9E67', Xmeso='#B090BF'),
               WT = c('early' = 'lightskyblue', 'mid' = 'steelblue3', 'late' = 'midnightblue'),
               Eed = c('1' = 'black', '0' = 'lightgrey'),
               Rnf2 = c('1' = 'black', '0' = 'lightgrey'),
               Kdm2b = c('1' = 'black', '0' = 'lightgrey')
    )
)

ho <- rowAnnotation(
    genes = anno_mark(at = which(!is.na(delta_expr$gene)),
    labels = delta_expr$gene[which(!is.na(delta_expr$gene))])
)

ht_expr <-
Heatmap(as.matrix(WT.g), show_row_names=F, cluster_rows=F, cluster_columns = diana(t(as.matrix(WT.g))), show_column_dend=F, col=colorRamp2(c(-2.5,0,2.5), c('magenta', 'black', 'yellow')), top_annotation = ha, right_annotation = ho)

ht_methyl + ht_expr + ht_diff

```


# ED 7c
```{r}
require(ggplot2)

data <- read.table('DMV_pos_percent.tsv', header=T)
data$time <- factor(data$time, levels=c('early','mid','late'))

ggplot(data, aes(x=time, fill=time, y=count)) + geom_col() + theme_classic() + scale_fill_manual(values=c('steelblue1', 'steelblue3', 'steelblue4'))

```













