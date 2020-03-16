# ED 8a
```{r}
require(ggplot2)

data <- read.table('Eed_Rnf2_CGImethylation_vs_WT.tsv', header=T)

ggplot(data, aes(x=Eed_Epi,y=Rnf2_Epi)) + theme_classic() + ylim(c(-0.5,0.5)) + xlim(c(-0.5,0.5)) + geom_vline(xintercept=c(-0.1,0,0.1)) + geom_hline(yintercept=c(-0.1,0,0.1)) + geom_point(size=1, aes(color=delta_Epi)) + scale_color_manual(values=c('cornflowerblue','lightsalmon','grey','cornflowerblue','blue','purple','lightsalmon','purple','red3')) + theme(legend.position='none') + theme(legend.position='none')

ggplot(data, aes(x=Eed_Exe,y=Rnf2_Exe)) + theme_classic() + ylim(c(-0.5,0.5)) + xlim(c(-0.5,0.5)) + geom_vline(xintercept=c(-0.1,0,0.1)) + geom_hline(yintercept=c(-0.1,0,0.1)) + geom_point(size=1, aes(color=delta_Exe)) + scale_color_manual(values=c('cornflowerblue','lightsalmon','grey','cornflowerblue','blue','purple','lightsalmon','purple','red3')) + theme(legend.position='none') + theme(legend.position='none')
```

