# Figure 3a
```{r}
require(ggplot2)

data <- read.table('Fig3a_PCG_count.tsv', header=T)

ggplot(data, aes(x=stage, y=frac)) + geom_boxplot() + theme_classic()
```


# Figure 3b
```{r}
require(ggplot2)

data <- read.table('Fig3b_Allantois_count.tsv', header=T)

ggplot(data, aes(x=stage, y=frac)) + geom_boxplot() + theme_classic()
```


# Figure 3d
```{r}
require(ggplot2)

data <- read.table('Fig3d_XEcto_fraction.tsv', header=T)

ggplot(data, aes(x=sex, y=frac)) + geom_point() + geom_boxplot() + theme_classic() + facet_grid(~type)
```


# Figure 3e
```{r}
require(ggplot2)

data <- read.table('Fig3e.ED6b.ED10f_X_expr.tsv', header=T)

ggplot(subset(data, stage %in% c('WT_65','WT_70','WT_75','WT_80','WT_85','Eed_85','Rnf2_85','Kdm2b_85')), aes(x=gsub('_.*','',stage), y=log2(frac), fill=sex)) + geom_boxplot(width=0.5, outlier.shape=NA) + theme_classic() + scale_fill_manual(values=c('firebrick','dodgerblue')) + facet_grid(~XE)
```


# Figure 3f
```{r}
require(ggplot2)

data <- read.table('Fig3f_Cdkn2a_pos_cells.tsv', header=T)

ggplot(data, aes(x=type, y=frac, fill=label)) + geom_boxplot(outlier.shape=NA) + theme_classic() + ylab('Fraction of Cdkn2a positive cells') + ylim(c(0,100)) + geom_point(position=position_jitterdodge(jitter.width=0.001), size=1)
```


# Figure 3h left
deeptools/bin/computeMatrix reference-point \
 -S Eed_Epi.bw Kdm2b_Epi_rep1.bw Rnf2_Epi_rep1.bw WT_Epi.bw \
 -R DMV_CGIs_cluster1.bed \
 --referencePoint center \
 -out DMV_CGIs_cluster1.tab.gz \
 -a 5000 -b 5000 \
 --binSize 100

deeptools/bin/plotProfile -m DMV_CGIs_cluster1.tab.gz -out DMV_CGIs_cluster1.pdf --averageType median --plotType std --yMin 0 --yMax 1

# Figure 3h right
deeptools/bin/computeMatrix reference-point \
 -S Eed_Exe.bw Kdm2b_Exe_rep1.bw Rnf2_Exe_rep1.bw WT_Exe.bw \
 -R DMV_CGIsExEhyper_cluster1.bed \
 --referencePoint center \
 -out DMV_CGIsExEhyper_cluster1.tab.gz \
 -a 5000 -b 5000 \
 --binSize 100 \
 -p 10

deeptools/bin/plotProfile -m DMV_CGIsExEhyper_cluster1.tab.gz -out DMV_CGIsExEhyper_cluster1.pdf --averageType median --plotType std --yMin 0 --yMax 1


# Figure 3i
```{r}
require(ggplot2)

WT_plot <- read.table('Fig3i_Dppa3_Umap_WT.tsv', header=T)
Eed_plot <- read.table('Fig3i_Dppa3_Umap_Eed.tsv', header=T)
Rnf2_plot <- read.table('Fig3i_Dppa3_Umap_Rnf2.tsv', header=T)
Kdm2b_plot <- read.table('Fig3i_Dppa3_Umap_Kdm2b.tsv', header=T)

ggplot() + geom_point(data=subset(WT_plot, X=='all'), aes(x=UMAP_1, y=UMAP_2), size=0.5, color='lightgrey') + geom_point(data=subset(WT_plot, X!='all'), aes(x=UMAP_1, y=UMAP_2, color=X), size=1) + theme_classic() + scale_color_manual(values=c('black','#E6007E')) + theme_void() + theme(legend.position="none") + ggtitle('')
ggplot() + geom_point(data=subset(Kdm2b_plot, X=='all'), aes(x=UMAP_1, y=UMAP_2), size=0.5, color='lightgrey') + geom_point(data=subset(Kdm2b_plot, X!='all'), aes(x=UMAP_1, y=UMAP_2, color=X), size=1) + theme_classic() + scale_color_manual(values=c('black','#E6007E')) + theme_void() + theme(legend.position="none") + ggtitle('')
ggplot() + geom_point(data=subset(Rnf2_plot, X=='all'), aes(x=UMAP_1, y=UMAP_2), size=0.5, color='lightgrey') + geom_point(data=subset(Rnf2_plot, X!='all'), aes(x=UMAP_1, y=UMAP_2, color=X), size=1) + theme_classic() + scale_color_manual(values=c('black','#E6007E')) + theme_void() + theme(legend.position="none") + ggtitle('')
ggplot() + geom_point(data=subset(Eed_plot, X=='all'), aes(x=UMAP_1, y=UMAP_2), size=0.5, color='lightgrey') + geom_point(data=subset(Eed_plot, X!='all'), aes(x=UMAP_1, y=UMAP_2, color=X), size=1) + theme_classic() + scale_color_manual(values=c('black','#E6007E')) + theme_void() + theme(legend.position="none") + ggtitle('')
```

