# ED 9b
```{r}
require(ggplot2)
require(EnrichedHeatmap)
require(circlize)

a <- 10

overview <- read.table('ED9b_Cutsites_E65_Eed_WT.tsv', header=T)
Eed <- read.table('ED9b_Cutsites_E65_Eed.tsv', header=T)
WT <- read.table('ED9b_Cutsites_E65_WT.tsv', header=T)

overview <- read.table('ED9b_Cutsites_E75_Eed_WT.tsv', header=T)
Eed <- read.table('ED9b_Cutsites_E75_Eed.tsv', header=T)
WT <- read.table('ED9b_Cutsites_E75_WT.tsv', header=T)

overview <- read.table('ED9b_Cutsites_E85_Eed_WT.tsv', header=T)
Eed <- read.table('ED9b_Cutsites_E85_Eed.tsv', header=T)
WT <- read.table('ED9b_Cutsites_E85_WT.tsv', header=T)

overview$Var1 <- factor(overview$Var1, levels=c('insertions','mismatch','match','spliced'))
ggplot(overview, aes(x=Var2, y=value, fill=Var1)) + geom_bar(stat='identity') + theme_classic() + facet_wrap(~type, ncol=1) + scale_fill_manual(values=c('darkblue','darkblue','gray90','darkgrey')) + xlab('Eed position [nt]') + ylab('Fraction of reads') + scale_x_continuous(breaks=c(cumsum(c(0,458,114,153,93,66,126,82,93,134,106,159,74,127,285)), 579,587,598, 939,950,961, 1052,1063,1074), labels=c(cumsum(c(0,458,114,153,93,66,126,82,93,134,106,159,74,127,285)), '','g3','', '','g2','', '','g1',''))


Eed_c <- Eed[,c((576-a):(598+a), (ncol(Eed)-1))]
WT_c <- WT[,c((576-a):(598+a), (ncol(WT)-1))]
Eed_c <- Eed_c[!apply(Eed_c[, -ncol(Eed_c)], 1, function(x){all(x=='.')}),]
WT_c <- WT_c[!apply(WT_c[, -ncol(WT_c)], 1, function(x){all(x=='.')}),]
ht_Eed_c <- Heatmap(Eed_c$EMBRYO, col=structure(2:20, names=levels(Eed_c$EMBRYO)), name='Embryos', show_row_names=FALSE, width=unit(3, 'mm')) +
Heatmap(as.matrix(Eed_c[,-ncol(Eed_c)]), show_row_names=F, show_column_names=T, cluster_rows=F, cluster_columns=F, col=structure(c('gray90','lightblue','lightblue','lightblue','lightblue','darkgrey','darkblue','darkblue','darkblue','darkblue','firebrick','orange'), names=c('.','A','T','G','C','-','a','t','g','c','I','n')), column_labels=gsub('\\..*','',colnames(Eed_c[,-ncol(Eed_c)])), column_names_rot=0)
ht_WT_c <- Heatmap(WT_c$EMBRYO, col=structure(2:20, names=levels(WT_c$EMBRYO)), name='Embryos', show_row_names=FALSE, width=unit(3, 'mm')) +
Heatmap(as.matrix(WT_c[,-ncol(WT_c)]), show_row_names=F, show_column_names=T, cluster_rows=F, cluster_columns=F, col=structure(c('gray90','lightblue','lightblue','lightblue','lightblue','darkgrey','darkblue','darkblue','darkblue','darkblue','firebrick','orange'), names=c('.','A','T','G','C','-','a','t','g','c','I','n')), column_labels=gsub('\\..*','',colnames(WT_c[,-ncol(WT_c)])), column_names_rot=0)

draw(ht_Eed_c, split=Eed_c$EMBRYO)
draw(ht_WT_c, split=WT_c$EMBRYO)



Eed_c <- Eed[,c((939-a):(961+a), (ncol(Eed)-1))]
WT_c <- WT[,c((939-a):(961+a), (ncol(Eed)-1))]
Eed_c <- Eed_c[!apply(Eed_c[, -ncol(Eed_c)], 1, function(x){all(x=='.')}),]
WT_c <- WT_c[!apply(WT_c[, -ncol(WT_c)], 1, function(x){all(x=='.')}),]
ht_Eed_c <- Heatmap(Eed_c$EMBRYO, col=structure(2:20, names=levels(Eed_c$EMBRYO)), name='Embryos', show_row_names=FALSE, width=unit(3, 'mm')) +
Heatmap(as.matrix(Eed_c[,-ncol(Eed_c)]), show_row_names=F, show_column_names=T, cluster_rows=F, cluster_columns=F, col=structure(c('gray90','lightblue','lightblue','lightblue','lightblue','darkgrey','darkblue','darkblue','darkblue','darkblue','firebrick','orange'), names=c('.','A','T','G','C','-','a','t','g','c','I','n')), column_labels=gsub('\\..*','',colnames(Eed_c[,-ncol(Eed_c)])), column_names_rot=0)
ht_WT_c <- Heatmap(WT_c$EMBRYO, col=structure(2:20, names=levels(WT_c$EMBRYO)), name='Embryos', show_row_names=FALSE, width=unit(3, 'mm')) +
Heatmap(as.matrix(WT_c[,-ncol(WT_c)]), show_row_names=F, show_column_names=T, cluster_rows=F, cluster_columns=F, col=structure(c('gray90','lightblue','lightblue','lightblue','lightblue','darkgrey','darkblue','darkblue','darkblue','darkblue','firebrick','orange'), names=c('.','A','T','G','C','-','a','t','g','c','I','n')), column_labels=gsub('\\..*','',colnames(WT_c[,-ncol(WT_c)])), column_names_rot=0)

draw(ht_Eed_c, split=Eed_c$EMBRYO)
draw(ht_WT_c, split=WT_c$EMBRYO)



Eed_c <- Eed[,c((1052-a):(1074+a), (ncol(Eed)-1))]
WT_c <- WT[,c((1052-a):(1074+a), (ncol(Eed)-1))]
Eed_c <- Eed_c[!apply(Eed_c[, -ncol(Eed_c)], 1, function(x){all(x=='.')}),]
WT_c <- WT_c[!apply(WT_c[, -ncol(WT_c)], 1, function(x){all(x=='.')}),]
ht_Eed_c <- Heatmap(Eed_c$EMBRYO, col=structure(2:20, names=levels(Eed_c$EMBRYO)), name='Embryos', show_row_names=FALSE, width=unit(3, 'mm')) +
Heatmap(as.matrix(Eed_c[,-ncol(Eed_c)]), show_row_names=F, show_column_names=T, cluster_rows=F, cluster_columns=F, col=structure(c('gray90','lightblue','lightblue','lightblue','lightblue','darkgrey','darkblue','darkblue','darkblue','darkblue','firebrick','orange'), names=c('.','A','T','G','C','-','a','t','g','c','I','n')), column_labels=gsub('\\..*','',colnames(Eed_c[,-ncol(Eed_c)])), column_names_rot=0)
ht_WT_c <- Heatmap(WT_c$EMBRYO, col=structure(2:20, names=levels(WT_c$EMBRYO)), name='Embryos', show_row_names=FALSE, width=unit(3, 'mm')) +
Heatmap(as.matrix(WT_c[,-ncol(WT_c)]), show_row_names=F, show_column_names=T, cluster_rows=F, cluster_columns=F, col=structure(c('gray90','lightblue','lightblue','lightblue','lightblue','darkgrey','darkblue','darkblue','darkblue','darkblue','firebrick','orange'), names=c('.','A','T','G','C','-','a','t','g','c','I','n')), column_labels=gsub('\\..*','',colnames(WT_c[,-ncol(WT_c)])), column_names_rot=0)

draw(ht_Eed_c, split=Eed_c$EMBRYO)
draw(ht_WT_c, split=WT_c$EMBRYO)
```

