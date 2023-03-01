#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(ggrepel)

## DEFINITIONS #####################################################################################

protein_quant_filename <- 'reports/protein_quantification.tsv'
ppmi_model_filename <- 'reports/PPMI2-1-model.tsv'
ppmi_data_filename <- 'reports/PPMI_data_before_model.tsv'

## FUNCTIONS #######################################################################################

ezwrite <- function(x, output_dir, output_filename) {
    # Wrapper for fwrite that uses standard TSV output defaults.
    # Concatenates output directory and filename for final output location.
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    fwrite(x, file=paste0(output_dir, '/', output_filename),
        quote=F,
        row.names=F,
        col.names=T,
        sep='\t')
}

## LOAD DATA #######################################################################################


ppmi.model <- fread(ppmi_model_filename)
ppmi.model <- ppmi.model[beta_value > 0 & rsquared >= 0.25]
ppmi.dat <- fread(ppmi_data_filename)




## Format for output for Hirotaka
if (FALSE) {
dat.csf <- merge(ppmi.dat[sampletype=='CSF'], ppmi.model[samples=='CSF'], by.x='UniProt', by.y='protein')
dat.pla <- merge(ppmi.dat[sampletype=='PLA'], ppmi.model[samples=='PLA'], by.x='UniProt', by.y='protein')
dat <- rbindlist(list(dat.csf, dat.pla))
dat <- dat[STUDY %in% c('PPMI1','PPMI2')]
dat[STUDY=='PPMI1', NPX := value]
dat[STUDY=='PPMI2', NPX := Intercept_value + (value*beta_value)]

dat[, 'ppea' := NULL]
dat[, 'rsquared_rank' := NULL]
dat[, 'Intercept_se' := NULL]
dat[, 'Intercept_t' := NULL]
dat[, 'Intercept_p' := NULL]
dat[, 'beta_se' := NULL]
dat[, 'beta_t' := NULL]
dat[, 'beta_p' := NULL]
dat[, 'sample' := NULL]
dat[STUDY=='PPMI1', 'Intercept_value' := NA]
dat[STUDY=='PPMI1', 'beta_value' := NA]
dat[STUDY=='PPMI1', 'rsquared' := NA]

setnames(dat, 'samples', 'SAMPLE')
setnames(dat, 'mode', 'olink_panel')
setnames(dat, 'NPX', 'NPX_modeled')
setnames(dat, 'value', 'NPX_original')
setnames(dat, 'samples', 'SAMPLE')

setcolorder(dat, c(
"UniProt",
"SAMPLE",
"olink_panel",
"ID",
"sampletype",
"month",
"BASELINE_AGE",
"SEX",
"PHENO",
"ARM",
"STUDY",
"Intercept_value",
"beta_value",
"rsquared",
"NPX_original",
"NPX_modeled"))

ezwrite(dat, 'reports/', 'modeled_NPX_values.tsv')
}



## APPLY SPECIFIC MODELS ###########################################################################
dat.csf <- merge(ppmi.dat[sampletype=='CSF'], ppmi.model[samples=='CSF'], by.x='UniProt', by.y='protein')
dat.pla <- merge(ppmi.dat[sampletype=='PLA'], ppmi.model[samples=='PLA'], by.x='UniProt', by.y='protein')
dat <- rbindlist(list(dat.csf, dat.pla))
dat <- dat[STUDY %in% c('PPMI1','PPMI2')]
dat[STUDY=='PPMI1', NPX := value]
dat[STUDY=='PPMI2', NPX := Intercept_value + (value*beta_value)]

g <- ggplot(dat, aes(x=value, y=NPX)) +
    geom_point(shape=21, alpha=0.3) +
    facet_grid(mode~sampletype) +
    theme_few() +
    labs(x='Original NPX', y='Model-adjusted NPX', title='CSF and PLA-specific model')

ggsave(g, file='figs/original-vs-model-NPX-specific.png', height=20, width=20, units='cm')

g2 <- ggplot(dat, aes(x=mode, y=NPX)) +
    geom_boxplot() +
    facet_grid(.~sampletype) +
    theme_few() +
    labs(x='Panel', y='Model-adjusted NPX', title='CSF and PLA-specific model') +
    geom_hline(yintercept=0, linetype='dashed', alpha=0.6, color='red') +
    coord_flip()

ggsave(g2, file='figs/modeled-NPX-distributions-specific.png', height=20, width=20, units='cm')

## CSF #############################################################################################
dat.csf <- merge(ppmi.dat[sampletype=='CSF'], ppmi.model[samples=='CSF_and_PLA'], by.x='UniProt', by.y='protein')
dat.pla <- merge(ppmi.dat[sampletype=='PLA'], ppmi.model[samples=='CSF_and_PLA'], by.x='UniProt', by.y='protein')
dat <- rbindlist(list(dat.csf, dat.pla))
dat <- dat[STUDY %in% c('PPMI1','PPMI2')]
dat[STUDY=='PPMI1', NPX := value]
dat[STUDY=='PPMI2', NPX := Intercept_value + (value*beta_value)]

g3 <- ggplot(dat, aes(x=value, y=NPX)) +
    geom_point(shape=21, alpha=0.3) +
    facet_grid(mode~sampletype) +
    theme_few() +
    labs(x='Original NPX', y='Model-adjusted NPX', title='Combined CSF+PLA model')

ggsave(g3, file='figs/original-vs-model-NPX-merged.png', height=20, width=20, units='cm')

g4 <- ggplot(dat, aes(x=mode, y=NPX)) +
    geom_boxplot() +
    facet_grid(.~sampletype) +
    theme_few() +
    labs(x='Panel', y='Model-adjusted NPX', title='Combined CSF+PLA model') +
    geom_hline(yintercept=0, linetype='dashed', alpha=0.6, color='red') +
    coord_flip()

ggsave(g4, file='figs/modeled-NPX-distributions-merged.png', height=20, width=20, units='cm')


## PLASMA ##########################################################################################

dat.pla