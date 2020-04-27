# Figure 4a
Staging (collected and corresponding adjusted) are in Supplementary Table 1.
```{r}
require(ggplot2)

data <- data.frame(collected=c(6.5, 6.5, 7.5, 8.5, 8.5), assigned=c(6.4, 7.5, 7.4, 7.0, 7.5), frac=c(15/17, 2/17, 7/7, 1/10, 9/10))

ggplot(data, aes(x=collected, y=0)) + geom_curve(aes(x=collected, y=0, xend=assigned, yend=0, size=frac), arrow=arrow(length = unit(0.01, 'npc')), curvature=1) + theme_classic()
```


# Figure 4c
Log2FC and expression correlation are in Fig4c_EED_FC_cor.tsv.


# Figure 4d
```{r}
require(eulerr)

data <- read.table('Fig4d.ED10g_Venn.tsv', header=T, sep='\t')

plot(euler(c(
    'WT_total'=2322,
    'WT_Dppa3'=8,
    'WT_Nanog'=975,
    'WT_Klf2'=63,
    'WT_Nanog&WT_Dppa3'=12,
    'WT_Klf2&WT_Dppa3'=0,
    'WT_Nanog&WT_Klf2'=71,
    'WT_Nanog&WT_Klf2&WT_Dppa3'=1,
    'EED_total'=2655,
    'EED_Dppa3'=47,
    'EED_Nanog'=972,
    'EED_Klf2'=163,
    'EED_Nanog&EED_Dppa3'=84,
    'EED_Klf2&EED_Dppa3'=15,
    'EED_Nanog&EED_Klf2'=188,
    'EED_Nanog&EED_Klf2&EED_Dppa3'=25
)), quantities=T, lables=list(font=4), shape='ellipse')
```

# Figure 4e
List of differentially expressed genes are in Supplementary Table 6.
Fraction positive cells for Hoxb1 and Hoxd9 are in Supplementary Table 9.


# Figure 4f
mESC Nanostring expression data are in Supplementary Table 11.
```{r}
require(ggplot2)
require(scales)

data <- read.table('Fig4f.ED11d_Nanostring.tsv', header=T)

data <- subset(data, gene %in% c('Dppa3', 'Esrrb', 'Pax6', 'Sox1', 'T', 'Tbx20', 'Bmp4', 'Gata4', 'Gata6'))
data$gene <- factor(data$gene, levels = rev(c('Dppa3', 'Esrrb', 'Pax6', 'Sox1', 'T', 'Tbx20', 'Bmp4', 'Gata4', 'Gata6')))
data$condition <- factor(data$condition, levels = c('2i','24hFgf', '72hFgf','RA', 'LDN', 'XAV', 'SB', 'B5', 'B500', 'W10', 'W100', 'W1000','A10', 'A1000','W.A'))

ggplot(data, aes(x=condition, y=gene)) +
  geom_tile(aes(fill = log2FC, color=sig)) +
  scale_fill_distiller(palette ='RdBu', direction = -1, na.value = 'grey', limits = c(-8,8), oob=squish) +
  scale_color_manual(values=c('grey', 'black'))
```

