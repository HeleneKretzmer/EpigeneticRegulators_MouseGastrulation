# ED 10b
```{r}
data <- read.table('ED10b_Embyro_sizes.tsv', header=T)

round(apply(data, 2, mean, na.rm=T),0)
round(apply(data, 2, sd, na.rm=T),0)
```


# ED 10c
```{r}
require(ggplot2)
require(ggalluvial)

cluster_colors <- c('black','#ffd100','#cbbba0','#92c464','#df76ac','#a2c5ea','#bd6fab','#477ec0','#5b57a2','#107f71','#00a19a','#3aaa35','#7a2182','#783f91','#f9b233','#e20146','#e94e1b','#9e937e','#b6a8d3','#ffed00','#581d6f','#6c6556','#30892b','#0070b2','#954b97','#5567ae','#b670ac','#8c1d82','#2fac66','#1c4024','#1a5b1a','#277424','#662483','#006f72','#6b3064','#66296b','#64358c','#009767','#642f2c','#4a712e','#f39200','#624758','#624758')
cluster_order <- c('17','4','25','8','37','3','2','38','36','11','26','19','13','20','29','27','14','0','5','7','9','28','1','30','12','23','15','34','35','33','24','39','22','31','41','18','32','16','21','10','40','6')
s <- c(17, 27, 8,19,1,11,16,10,24,39,33,35,26,31, 3,7,38,23, 37,2,13,20,12,30,36,32,34,21,6,18,9,22, 5,15,41, 0,25,28, 4,29,40,14)

data.WT <- read.table('ED10c_WT_alluvial.tsv', header=T)
data.Eed <- read.table('ED10c_Eed_alluvial.tsv', header=T)

data.WT$cluster <- factor(data.WT$cluster, levels=cluster_order[match(s, cluster_order)])
data.Eed$cluster <- factor(data.Eed$cluster, levels=cluster_order[match(s, cluster_order)])

ggplot(data.WT, aes(x = stage, stratum = cluster, alluvium = cluster, y = fraction, fill = cluster, label = cluster)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(alpha = .5) +
  geom_stratum(color = 'lightgrey') +
  theme_classic() +
  scale_fill_manual(values=cluster_colors[match(s, cluster_order)])

ggplot(data.Eed, aes(x = stage, stratum = cluster, alluvium = cluster, y = fraction, fill = cluster, label = cluster)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(alpha = .5) +
  geom_stratum(color = 'lightgrey') +
  theme_classic() +
  scale_fill_manual(values=cluster_colors[match(s, cluster_order)])
```


# ED 10d
```{r}
require(ggplot2)

data <- read.table('ED10d_Estimated_number_of_PGC.tsv', header=T)

ggplot(data, aes(x=adjStage, y=PGCs.estimate)) + geom_point() + geom_boxplot(fill=NA) + theme_classic() + facet_grid(~type)
```


# ED 10e
```{r}
require(ggplot2)

data <- read.table('ED10e_Cdkn2a_pos_cells_Eed.tsv', header=T)

ggplot(data, aes(x=type, y=frac, fill=label)) + geom_boxplot(outlier.shape=NA) + theme_classic() + ylab('Fraction of Cdkn2a positive cells') + ylim(c(0,100)) + geom_point(position=position_jitterdodge(jitter.width=0.001), size=1)
```


# ED 10f
```{r}
require(ggplot2)

data <- read.table('Fig3e.ED6b.ED10f_X_expr.tsv', header=T)

ggplot(subset(data, stage %in% c('WT_65','WT_70','WT_75','WT_80','WT_85','Eed_65','Eed_75','Eed_85')), aes(x=stage, y=log2(frac), fill=sex)) + geom_boxplot(width=0.5, outlier.shape=NA) + theme_classic() + scale_fill_manual(values=c('firebrick','dodgerblue')) + ggtitle('EED') + facet_grid(~XE)
```


# ED 10g
```{r}
require(eulerr)

data <- read.table('Fig4d.ED10g_Venn.tsv", header=T, sep='\t')

plot(euler(c(
    'WT_total'=4808,
    'WT_Dppa3'=42,
    'WT_Nanog'=1100,
    'WT_Tbx3'=67,
    'WT_Nanog&WT_Dppa3'=23,
    'WT_Tbx3&WT_Dppa3'=0,
    'WT_Nanog&WT_Tbx3'=20,
    'WT_Nanog&WT_Tbx3&WT_Dppa3'=1,
    'EED_total'=5335,
    'EED_Dppa3'=511,
    'EED_Nanog'=1462,
    'EED_Tbx3'=86,
    'EED_Nanog&EED_Dppa3'=365,
    'EED_Tbx3&EED_Dppa3'=13,
    'EED_Nanog&EED_Tbx3'=68,
    'EED_Nanog&EED_Tbx3&EED_Dppa3'=16
)), quantities=T, lables=list(font=4), shape='ellipse')

plot(euler(c(
    'WT_total'=787,
    'WT_Dppa3'=31,
    'WT_Nanog'=295,
    'WT_Klf2'=18,
    'WT_Nanog&WT_Dppa3'=26,
    'WT_Klf2&WT_Dppa3'=0,
    'WT_Nanog&WT_Klf2'=8,
    'WT_Nanog&WT_Klf2&WT_Dppa3'=5,
    'EED_total'=2044,
    'EED_Dppa3'=150,
    'EED_Nanog'=517,
    'EED_Klf2'=162,
    'EED_Nanog&EED_Dppa3'=179,
    'EED_Klf2&EED_Dppa3'=76,
    'EED_Nanog&EED_Klf2'=187,
    'EED_Nanog&EED_Klf2&EED_Dppa3'=96
)), quantities=T, lables=list(font=4), shape='ellipse')
```















