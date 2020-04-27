# ED 11c
```{r}
require(ggplot2)
require(scales)

data <- read.table('Fig4f.ED11c_Nanostring.tsv', header=T)

data$gene <- factor(data$gene, levels = rev(c('Dppa3', 'Esrrb', 'Prdm1', 'Prdm14', 'Pou5F1', 'Nanos3', 'Sox2', 'Nanog', 'Dnd1', 'Pax6', 'Sox1', 'Nes', 'Dlx5', 'Eya', 'Sox3', 'Foxd3', 'Eomes', 'T', 'Mixl1', 'Evx1', 'Lhx1', 'Mesp1', 'Dll1', 'Tbx6', 'Noto', 'Shh', 'Gsc', 'Kdr', 'Foxf1', 'Tbx4', 'Tbx20', 'Postn', 'Gata4', 'Gata6', 'Foxa2', 'Sox17', 'Elf5', 'Bmp4', 'Bmp8b', 'Fgf5', 'Fgf4', 'Wnt3a', 'Nodal', 'Notch1')))
data$condition <- factor(data$condition, levels = c('2i','24hFgf', '72hFgf','RA', 'LDN', 'XAV', 'SB', 'B5', 'B500', 'W10', 'W100', 'W1000','A10', 'A1000','W.A'))

ggplot(data, aes(x=condition, y=gene)) +
  geom_tile(aes(fill = log2meanexpr.x, color=sig)) +
  scale_fill_distiller(palette ='RdBu', direction = -1, na.value = 'grey', limits = c(0,12), oob=squish) +
  scale_color_manual(values=c('grey', 'black')) + 
  ggtitle('V65 (WT mESC line)')

ggplot(data, aes(x=condition, y=gene)) +
  geom_tile(aes(fill = log2meanexpr.y, color=sig)) +
  scale_fill_distiller(palette ='RdBu', direction = -1, na.value = 'grey', limits = c(0,12), oob=squish) +
  scale_color_manual(values=c('grey', 'black')) + 
  ggtitle('c22 (EED KO mESC line)')
```


# ED 11d
```{r}
require(ggplot2)
require(scales)
require(ggsignif)

data <- read.table('ED11d_Nanostring_means.tsv', header=T)
sig <- read.table('ED11d_Nanostring_sig.tsv', header=T)

data$gene <- factor(data$gene, levels = c('Dppa3','Esrrb','Pax6', 'Sox1','T','Tbx20','Bmp4', 'Bmp8b','Gata4', 'Gata6'))
data$condition <- factor(data$condition, levels = c('2i','24hFgf', '72hFgf','RA', 'LDN', 'XAV', 'SB', 'B5', 'B500', 'W10', 'W100', 'W1000','A10', 'A1000','W.A'))

sig$gene <- factor(sig$gene, levels = c('Dppa3','Esrrb','Pax6', 'Sox1','T','Tbx20','Bmp4', 'Bmp8b','Gata4', 'Gata6'))
sig$condition <- factor(sig$condition, levels = c('2i','24hFgf', '72hFgf','RA', 'LDN', 'XAV', 'SB', 'B5', 'B500', 'W10', 'W100', 'W1000','A10', 'A1000','W.A'))

ggplot(data, aes(x=type, y=mean_expr, fill=type)) + 
  geom_col(position='dodge') +
  geom_errorbar(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr), width=0.2) + 
  facet_grid(gene~condition, scales='free_y') +
  theme_classic() + scale_y_continuous(expand = c(0,0,0.25,0)) +
  geom_signif(data = sig, aes(x=type, y=expression), comparisons = list(c('V65', 'c22')), test='t.test', map_signif_level=F, textsize = 2)
```

