# ED 4b
```{r}
require(pheatmap)

data <- read.table('ED4a.ED4b_WGBS_genomewide.tsv', header=T)

M <- cor(data[,-c(1:3)], use='pairwise.complete.obs')

pheatmap(M[grep('Epi', rownames(M)),grep('Epi', colnames(M))])
pheatmap(M[grep('Ex', rownames(M)),grep('Ex', colnames(M))])
```



# ED 4c
```{r}
require(vioplot)

data <- read.table('ED4a.ED4b_WGBS_genomewide.tsv', header=T)

vioplot(data[,-c(1:3)], names=c(colnames(data[,-c(1:3)])), cex=0.3)
```


# ED 4d
```{r}
require(ggplot)

data <- read.table('ED4c_CGI_avg_methyl.tsv', header=T)

ggplot(data, aes(x=WT_Epi, y=Eed_Epi)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Eed_Epi-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=Eed_Exe)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Eed_Exe-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')

ggplot(data, aes(x=WT_Epi, y=Rnf2_Epi_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Rnf2_Epi_rep1-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=Rnf2_Exe_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Rnf2_Exe_rep1-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')

ggplot(data, aes(x=WT_Epi, y=Kdm2b_Epi_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Kdm2b_Epi_rep1-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=Kdm2b_Exe_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Kdm2b_Exe_rep1-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')

ggplot(data, aes(x=WT_Epi, y=L3Mbtl2_Epi_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(L3Mbtl2_Epi_rep1-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=L3Mbtl2_Exe_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(L3Mbtl2_Exe_rep1-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')

ggplot(data, aes(x=WT_Epi, y=Dnmt1_Epi_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Dnmt1_Epi_rep1-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=Dnmt1_Exe_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Dnmt1_Exe_rep1-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')

ggplot(data, aes(x=WT_Epi, y=Dnmt3a_Epi_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Dnmt3a_Epi_rep1-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=Dnmt3a_Exe_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Dnmt3a_Exe_rep1-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')

ggplot(data, aes(x=WT_Epi, y=Dnmt3b_Epi_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Dnmt3b_Epi_rep1-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=Dnmt3b_Exe_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Dnmt3b_Exe_rep1-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')

ggplot(data, aes(x=WT_Epi, y=G9a_Epi_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(G9a_Epi_rep1-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=G9a_Exe_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(G9a_Exe_rep1-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')

ggplot(data, aes(x=WT_Epi, y=Kmt2a_Epi_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Kmt2a_Epi_rep1-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=Kmt2a_Exe_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Kmt2a_Exe_rep1-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')

ggplot(data, aes(x=WT_Epi, y=Kmt2b_Epi_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Kmt2b_Epi_rep1-WT_Epi, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
ggplot(data, aes(x=WT_Exe, y=Kmt2b_Exe_rep1)) + geom_abline(intercept=0, slope=1, color='grey') + geom_point(aes(color=cut(Kmt2b_Exe_rep1-WT_Exe, breaks=c(-1,-0.25,-0.1,0.1,0.25,1)))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position='none')
```


# ED 4e left
```{r}
require(eulerr)

data <- read.table('ED4e_Promoter_CGI_Kdm2b_Kmt2b.tsv', header=T)
data <- na.omit(data)
data <- data[rowSums(data[,-1]>0.1)>0,]
data$Kdm2b_Epi_up <- ifelse(data$Kdm2b_Epi_delta > 0.1, 'Kdm2b_Epi_up', 'Kdm2b_Epi_not')
data$Kmt2b_Epi_up <- ifelse(data$Kmt2b_Epi_delta > 0.1, 'Kmt2b_Epi_up', 'Kmt2b_Epi_not')
data$Kdm2b_Exe_up <- ifelse(data$Kdm2b_Exe_delta > 0.1, 'Kdm2b_Exe_up', 'Kdm2b_Exe_not')
data$Kmt2b_Exe_up <- ifelse(data$Kmt2b_Exe_delta > 0.1, 'Kmt2b_Exe_up', 'Kmt2b_Exe_not')

table(apply(data[,c(6:9)], 1, paste, collapse=','))

plot(euler(c(
    'Kmt2b_Exe'=781,
    'Kdm2b_Exe'=4,
    'Kdm2b_Exe&Kmt2b_Exe'=3,
    'Kmt2b_Epi'=384,
    'Kmt2b_Epi&Kmt2b_Exe'=390,
    'Kmt2b_Epi&Kdm2b_Exe'=1,
    'Kmt2b_Epi&Kdm2b_Exe&Kmt2b_Exe'=2,
    'Kdm2b_Epi'=55,
    'Kdm2b_Epi&Kmt2b_Exe'=13,
    'Kdm2b_Epi&Kmt2b_Epi'=99,
    'Kdm2b_Epi&Kmt2b_Epi&Kmt2b_Exe'=77,
    'Kdm2b_Epi&Kmt2b_Epi&Kdm2b_Exe'=1,
    'Kdm2b_Epi&Kmt2b_Epi&Kdm2b_Exe&Kmt2b_Exe'=4
)), quantities=T, lables=list(font=4), shape='ellipse')
```


# ED 4e right
```{r}
require(ggplot2)

data <- read.table('ED4e_Promoter_CGI_Kdm2b_Kmt2b_expr.tsv', header=T)

ggplot(data, aes(x=stage, y=perc.pos.cells)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=label)) + theme_classic() + ylab('Fraction of positive cells') + ylim(c(0,100)) + facet_grid(~delta)
```


# ED 4f
```{r}
require(ggplot2)

data <- read.table('ED4f_Distance_diff.methyl.CpG_to_CGI.tsv', header=T)

ggplot(data, aes(x=dist, y=frac, color=KO, group=KO)) + geom_line(size=1) + theme_classic() + xlim(c(0,10000)) + ylim(c(0,0.025)) + scale_color_manual(values=c('black','grey90','grey60','red','grey30')) + ylab('Distribution of hypermethylated CpGs') + xlab('Distance to CGI center [nt]')
```



