#!/usr/bin/env Rscript
# install.packages('OlinkAnalyze', repos='http://cran.us.r-project.org')"

# library(OlinkAnalyze)
library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)

protein_quant_filename <- 'reports/protein_quantification.tsv'
ppmi_model_filename <- 'reports/PPMI2-1-model.tsv'
load_protein_data <- function() {
    if(file.exists(protein_quant_filename)) {
        dat <- fread(protein_quant_filename)
        return(dat)
    } else {

        dd <- '/data/CARD/projects/proteomicsXprogression_PD/data-allreleases/amp-pd-v25-proteomics/'
        filenames <- list.files(path=dd, pattern="NPX.*csv$", recursive=T, full.names=T)

        # combine data sets and transform to long format
        dat <- foreach(filename=filenames, .combine='rbind') %do% {
                dat.tmp <- fread(filename)
                tmpname <- strsplit(strsplit(basename(filename), split='PreviewRelease_|_NPX\\.csv')[[1]][2], split='_')[[1]]
                mode <- tmpname[1]
                sampletype <- tmpname[2]
                dat.tmp.long <- melt(dat.tmp, 
                                measure.vars=colnames(dat.tmp)[2:ncol(dat.tmp)],
                                variable.name='sample')
                dat.tmp.long[, 'mode' := mode]
                return(dat.tmp.long[])
                }

        # split 'variable' column into its separate values
        dat[, c('PD_or_PP', 'participant_id', 'visit', 'sampletype', 'ppea', 'D') :=
            tstrsplit(sample, split='-')]

        # remove redundant columns
        dat[, sample := NULL]
        dat[, ppea := NULL]


        # Simplify visit codes to numeric
        dat[visit=="BLM0T1",  month:= 0 ]
        dat[visit=="SVM3T1",  month:= 3 ]
        dat[visit=="SVM6T1",  month:= 6 ]
        dat[visit=="SVM9T1",  month:= 9 ]
        dat[visit=="SVM12T1", month:= 12]
        dat[visit=="SVM18T1", month:= 18]
        dat[visit=="SVM24T1", month:= 24]
        dat[visit=="SVM30T1", month:= 30]
        dat[visit=="SVM36T1", month:= 36]
        dat[visit=="SVM42T1", month:= 42]
        dat[visit=="SVM48T1", month:= 48]
        dat[visit=="SVM54T1", month:= 54]
        dat[visit=="SVM60T1", month:= 60]
        dat[visit=="SVM72T1", month:= 72]
        dat[visit=="SVM84T1", month:= 84]
        dat[visit=="SVM96T1", month:= 96]
        dat[, visit := NULL]

        # Simplify D02 to 2, D01 to 1
        dat[D=='D02', D := '2']
        dat[D=='D01', D := '1']

        fwrite(dat, file=protein_quant_filename, quote=F, row.names=F, col.names=T, sep='\t')
        return(dat)
    }
}


load_enrollment_data <- function(cohort) {
    # load enrollment data (phenotypes at enrollment)
    enrollment_fn <- paste0('/data/CARD/projects/proteomicsXprogression_PD/data-allreleases/',
                            'amppd-v25-baseline-casecon-age-sex-enrollment.txt')
    enrollment <- fread(enrollment_fn, header=TRUE)
    enrollment <- rbindlist(list(
                    enrollment[PHENO==0 & ENROLL_STUDY_ARM=='LBD'],
                    enrollment[PHENO==1 & ENROLL_STUDY_ARM=='Healthy Control'],
                    enrollment[PHENO==2 & ENROLL_STUDY_ARM=='PD']
                ))
    dat.out <- enrollment[COHORT == cohort][]
    setkey(dat.out, ID)
    return(dat.out)
}

## LOAD AND FORMAT #################################################################################
dat <- load_protein_data()
dat[, ID := paste0(PD_or_PP, '-', participant_id)]
setkey(dat, ID)

ppmi_ids <- load_enrollment_data('PPMI')
ppmi <- merge(dat, ppmi_ids)

pdbp_ids <- load_enrollment_data('PDBP')
pdbp <- merge(dat, pdbp_ids)

dat.all <- rbindlist(list(ppmi, pdbp))
dat.all[participant_id == 'PDZZ318BWK', .N, by=COHORT]
dat.all[, STUDY := paste0(COHORT, D)]

dat.all[, D := NULL]
dat.all[, participant_id := NULL]
dat.all[, COHORT := NULL]
dat.all[, PD_or_PP := NULL]
setnames(dat.all, 'ENROLL_STUDY_ARM', 'ARM')

PPMI1_ids <- dat.all[STUDY=='PPMI1', ID]
PPMI2_ids <- dat.all[STUDY=='PPMI2', ID]

in_both_ppmi <- intersect(PPMI1_ids, PPMI2_ids)

dat.bothppmi <- dat.all[ID %in% in_both_ppmi]

# Reshape to wide to have paired observations, looking only at month 0
dat.bothppmi.wide <- dcast(dat.bothppmi[month==0], ID+UniProt+mode+sampletype+month+ARM~STUDY, value.var='value')
dat.bothppmi.wide <- dat.bothppmi.wide[! is.na(PPMI1)][! is.na(PPMI2)]



# Modeling PPMI D02 to D01
build_model <- function(DT, UniProt.i, sampletype.i) {
    # Subset to protein name and sample type of interest
    DT.sub <- DT[UniProt == UniProt.i][sampletype %in% sampletype.i]

    # Run linear model
    mod <- lm(data=DT.sub, PPMI1 ~ PPMI2)

    # Build output table from summary stats
    res <- data.table(coef(summary(mod)), keep.rownames=T)
    setnames(res, c('V','value','se','t','p'))                  # Simplify output col names
    res[, 'protein' := UniProt.i]
    res[V=='(Intercept)', V := 'Intercept']                     # Rename without parantheses
    res[V=='PPMI2', V := 'beta']                                # Rename to 'beta'

    # Reshape to long (in prep for wide format)
    res.long <- melt(res, measure.vars=c('value','se','t','p'))

    # Reshape to wide (single row)
    out <- dcast(res.long, protein~V+variable, value.name='value')

    # add r-squared column
    out[, 'rsquared' := summary(mod)$r.squared]

    return(out)
}


# Generate model summary, or load if it already exists
if(! file.exists(ppmi_model_filename)) {
    # summarize for each sample combo
    o <- foreach(sampletype.i=list('CSF','PLA',c('CSF','PLA')), .combine='rbind', .errorhandling='pass') %do% {
        foreach(UniProt.i=unique(dat.test$UniProt), .combine='rbind', .errorhandling='pass') %do% {

                out <- build_model(dat.bothppmi.wide, UniProt.i, sampletype.i)
                out[, samples := paste(sampletype.i, collapse='_and_')]
                return(out)

        }
    }
    o[, rsquared_rank := frank(-rsquared), by=samples]
    fwrite(o, file=ppmi_model_filename, quote=F, row.names=F, col.names=T, sep='\t')
} else {
    o <- fread(ppmi_model_filename)
}

# Convert to long format for plotting
o.long <- melt(o, measure.vars=c('Intercept_value','beta_value', 'rsquared'))
o.long[, independent_rank := frank(-value), by=list(variable, samples)]


# ggplot labelle function
labeller <- c(
    `CSF`='model: CSF Only',
    `PLA`='model: PLA Only',
    `CSF_and_PLA`='model: Combined CSF and PLA',
    `Intercept_value`='Intercept',
    `beta_value`='beta',
    `rsquared`='R-squared'
)


# Plot ranked by r-squared alone
g.ranks.rsquared <-  ggplot(o.long, aes(x=rsquared_rank, y=value, color=samples)) + 
    geom_point(shape=21, alpha=0.5) +
    facet_grid(variable~., scales='free', labeller=as_labeller(labeller), switch='both') +
    scale_color_discrete(name='Samples used in model') +
    labs(x='Rank (all facets ranked by descending R-squared to show\nhow beta/intercept change with R-squared)')

# rank independently for intercept, beta, rsquared)
g.ranks.independent <- ggplot(o.long, aes(x=independent_rank, y=value, color=samples)) +
    geom_point(shape=21, alpha=0.5) +
    scale_color_discrete(name='Samples used in model') +
    facet_grid(variable~., scales='free', labeller=as_labeller(labeller), switch='both') +
    labs(x='Ranked (all facets ranked separately, to show distribution')

# rank independently for intercept, beta, rsquared)
g.ranks.beeswarm <- ggplot(o.long, aes(x=1, y=value)) +
    geom_beeswarm() +
    facet_grid(variable~samples, scales='free', labeller=as_labeller(labeller), switch='both') +
    scale_x_continuous(breaks=NULL, labels=NULL) +
    labs(x=NULL)


ggsave(g.ranks.rsquared, file='figs/ppmi-ranks-rsquared.jpg', width=25, height=25, units='cm')
ggsave(g.ranks.independent, file='figs/ppmi-ranks-independent.jpg', width=25, height=25, units='cm')
ggsave(g.ranks.beeswarm, file='figs/ppmi-ranks-beeswarm.jpg', width=25, height=25, units='cm')


quit()