#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(ggrepel)

files <- list.files(pattern='.tsv')

interactome <- unique(c('CEP55','PDCD6','TFG','PDCD6','WWOX','FUBP1',
                    'CEP55','TFG','NFKBID','KRTAP13-2','ARSA',
                    'KRTAP6-2','PLSCR1','P36810	','CST6','P08563-PRO_0000041303',
                    'ENO2','HNRNPH3','A0A380PL96','A0A6L8PQC3','ZDHHC17','Q5NHS4',
                    'CDC5L','EGFR','PRKN','CD81','CUL5','JUN','CD81','ZDHHC5','HTRA4','TP53BP1'))

o2 <- foreach(file=files, .combine='rbind') %do% {
    dat <- fread(file)
    fn <- strsplit(file, split='ANXA11|_|-|\\.')[[1]]
    t1 <- fn[3]
    t2 <- fn[7]
    treatment <- paste0(t1, '_vs_', t2)

    o <- foreach(i=interactome, .combine='rbind') %do% {
        dat[Genes %like% i]
        }

    o <- unique(o)
    o[, 'treatment' := treatment]

}


o2[p.adj > 0.05, Genes := NA]


g <- ggplot(o2, aes(label=Genes, x=log2(ratio), y=-1*log10(p.adj))) + 
geom_label_repel() + 
geom_point() + 
facet_wrap(~treatment) +
geom_hline(yintercept=1.3, linetype='dashed', alpha=0.5)

ggsave(g, file='interactome.png', width=25, height=25, units='cm')