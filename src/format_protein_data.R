#!/usr/bin/env Rscript
# install.packages('OlinkAnalyze', repos='http://cran.us.r-project.org')"

library(OlinkAnalyze)
library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)

protein_quant_filename <- 'reports/protein_quantification.tsv'

load_protein_data <- function() {
    if(file.exists(protein_quant_filename)) {
        dat <- fread(protein_quant_filename)
        return(dat)
    } else {

        dd <- '/data/CARD/projects/proteomicsXprogression_PD/data-allreleases/amp-pd-v3-proteomics/'
        filenames <- list.files(path=dd, pattern="D0[12]_matrix_.*csv$", recursive=T, full.names=T)

        # combine data sets and transform to long format
        dat <- foreach(filename=filenames, .combine='rbind') %do% {
                mode <- strsplit(basename(filename), split='_matrix_|\\.csv')[[1]][2]
                dat.tmp <- fread(filename)
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


## OK ##############################################################################################
dat.all[, 'ppea' := NULL]
ppmi1 <- unique(dat.all[STUDY=='PPMI1']$ID)
ppmi2 <- unique(dat.all[STUDY=='PPMI2']$ID)

bridges <- intersect(ppmi1, ppmi2)

ppmi1.dat <- dat.all[STUDY=='PPMI1']
ppmi2.dat <- dat.all[STUDY=='PPMI2']

setkey(ppmi1.dat, ID, sampletype, month, mode, UniProt, ARM, PHENO, SEX)
setkey(ppmi2.dat, ID, sampletype, month, mode, UniProt, ARM, PHENO, SEX)
ppmi.merge <- merge(ppmi1.dat, ppmi2.dat)

ggplot(ppmi.merge[PHENO==1], aes(x=value.x, y=value.y, color=mode, shape=sampletype)) + geom_point()


in_both_ppmi <- intersect(ppmi1, ppmi2)

# Check M4 
M4 <- fread('reports/M4.tsv')

## 39/186 of the M4 genes are containd within the olink panel
#  UniProt           panel
#   P00568 cardiometabolic
#   P04080 cardiometabolic
#   P04792 cardiometabolic
#   P05556 cardiometabolic
#   P07339 cardiometabolic
#   P17931 cardiometabolic
#   Q04760 cardiometabolic
#   Q9UBP4 cardiometabolic
#   P36959    inflammation
#   P52564    inflammation
#   Q12765    inflammation
#   Q9H008    inflammation
#   Q9Y2J8    inflammation
#   O95817       neurology
#   P06733       neurology
#   P07108       neurology
#   P08758       neurology
#   P09211       neurology
#   P14868       neurology
#   P15311       neurology
#   P29218       neurology
#   P30043       neurology
#   P50135       neurology
#   Q01469       neurology
#   Q06323       neurology
#   Q06830       neurology
#   Q96IU4       neurology
#   Q99497       neurology
#   O14558        oncology
#   O94760        oncology
#   P09960        oncology
#   P14136        oncology
#   P15121        oncology
#   P30041        oncology
#   P40121        oncology
#   P43490        oncology
#   P51858        oncology
#   Q00796        oncology
#   Q99536        oncology

# Merge in enrollment status

# Look at delta-M4 vs enrollment status
dat.all[, 'M4' := ifelse(UniProt %in% M4$UniProt, TRUE, FALSE)]

# look at data that exists in quadruplicate
quad_samples <- copy(dat.all[, .N, by=list(ID, UniProt, sampletype, month, STUDY)][N ==4, !c('N')])
setkey(quad_samples, ID, UniProt, sampletype, month, STUDY)
setkey(dat.all, ID, UniProt, sampletype, month, STUDY)
quad_samples <- merge(dat.all, quad_samples)
quad.wide <- dcast(quad_samples, ID+UniProt+sampletype+month+STUDY~mode, value.var='value')
quad.sub <- quad.wide[month %in% c(0,12,24,36)]


if(FALSE) {
c_v_i <- ggplot(data=quad.sub) +
    geom_abline(slope=1, intercept=0, linetype='dashed', alpha=0.3) +
    geom_point(alpha=0.3, aes(x=cardiometabolic, y=inflammation, color=STUDY, shape=UniProt)) +
    facet_grid(sampletype~month) +
    guides(color = guide_legend(override.aes = list(shape=15, alpha=1, size=3))) +
    theme_few(12) +
    scale_x_continuous(breaks=0, labels=NULL) +
    scale_y_continuous(breaks=0, labels=NULL)

c_v_n <- ggplot(data=quad.sub) +
    geom_abline(slope=1, intercept=0, linetype='dashed', alpha=0.3) +
    geom_point(alpha=0.3, aes(x=cardiometabolic, y=neurology, color=STUDY, shape=UniProt)) +
    facet_grid(sampletype~month) +
    guides(color = guide_legend(override.aes = list(shape=15, alpha=1, size=3))) +
    theme_few(12) +
    scale_x_continuous(breaks=0, labels=NULL) +
    scale_y_continuous(breaks=0, labels=NULL)

c_v_o <- ggplot(data=quad.sub) +
    geom_abline(slope=1, intercept=0, linetype='dashed', alpha=0.3) +
    geom_point(alpha=0.3, aes(x=cardiometabolic, y=oncology, color=STUDY, shape=UniProt)) +
    facet_grid(sampletype~month) +
    guides(color = guide_legend(override.aes = list(shape=15, alpha=1, size=3))) +
    theme_few(12) +
    scale_x_continuous(breaks=0, labels=NULL) +
    scale_y_continuous(breaks=0, labels=NULL)

i_v_n <- ggplot(data=quad.sub) +
    geom_abline(slope=1, intercept=0, linetype='dashed', alpha=0.3) +
    geom_point(alpha=0.3, aes(x=inflammation, y=neurology, color=STUDY, shape=UniProt)) +
    facet_grid(sampletype~month) +
    guides(color = guide_legend(override.aes = list(shape=15, alpha=1, size=3))) +
    theme_few(12) +
    scale_x_continuous(breaks=0, labels=NULL) +
    scale_y_continuous(breaks=0, labels=NULL)

i_v_o <- ggplot(data=quad.sub) +
    geom_abline(slope=1, intercept=0, linetype='dashed', alpha=0.3) +
    geom_point(alpha=0.3, aes(x=inflammation, y=oncology, color=STUDY, shape=UniProt)) +
    facet_grid(sampletype~month) +
    guides(color = guide_legend(override.aes = list(shape=15, alpha=1, size=3))) +
    theme_few(12) +
    scale_x_continuous(breaks=0, labels=NULL) +
    scale_y_continuous(breaks=0, labels=NULL)

n_v_o <- ggplot(data=quad.sub) +
    geom_abline(slope=1, intercept=0, linetype='dashed', alpha=0.3) +
    geom_point(alpha=0.3, aes(x=neurology, y=oncology, color=STUDY, shape=UniProt)) +
    facet_grid(sampletype~month) +
    guides(color = guide_legend(override.aes = list(shape=15, alpha=1, size=3))) +
    theme_few(12) +
    scale_x_continuous(breaks=0, labels=NULL) +
    scale_y_continuous(breaks=0, labels=NULL)

ggsave(c_v_i, file='figs/cardio-v-inflammation.png', width=20, height=10, units='cm')
ggsave(c_v_n, file='figs/cardio-v-neuro.png', width=20, height=10, units='cm')
ggsave(c_v_o, file='figs/cardio-v-onco.png', width=20, height=10, units='cm')
ggsave(i_v_n, file='figs/inflammation-v-neuro.png', width=20, height=10, units='cm')
ggsave(i_v_o, file='figs/inflammation-v-onco.png', width=20, height=10, units='cm')
ggsave(n_v_o, file='figs/neuro-v-onco.png', width=20, height=10, units='cm')
}


# randomly select one value for duplicated proteins
set.seed(100)
dat.all <- dat.all[dat.all[, .I[sample(.N, size=1)], by=list(ID, UniProt, sampletype, month, STUDY)]$V1]
dat.wide <- dcast(dat.all, ARM+M4+ID+UniProt+month+STUDY~sampletype, value.var='value')

# Look at M4 proteins
labeller <- c(
    `TRUE`='After Baseline',
    `FALSE`='Baseline',
    `Healthy Control`='Healthy Control',
    `PD`='PD Diagnosis'
)
g.M4_A <- ggplot(dat.wide[month %in% c(0,12,24,36)], aes(x=CSF, y=PLA, color=M4)) +
    geom_point(shape=21, alpha=0.1) +
    facet_grid(month>0~ARM, labeller=as_labeller(labeller)) +
    scale_color_discrete(name='M4 Module') +
    labs(x='CSF', y='Plasma') +
    guides(color = guide_legend(override.aes=c(shape=16, size=3, alpha=1)))
ggsave(g.M4_A, file='figs/M4_A.jpg', width=20, height=20, units='cm')


labeller <- c(
    `TRUE`='in M4 module',
    `FALSE`='not in M4 module',
    `Healthy Control`='Healthy Control',
    `PD`='PD Diagnosis'
)
dat.all[, month := as.factor(month)]
g.M4_B <- ggplot(dat.all[month %in% c(0,24)], aes(x=month, y=value, color=sampletype)) +
    geom_boxplot() +
    scale_color_discrete(
            name='Sample Type',
            labels=c(`CSF`='CSF', `PLA`='Plasma')) +
    labs(x='Month', y='Olink Quantification Value') +
    facet_grid(M4~ARM, labeller=as_labeller(labeller))
ggsave(g.M4_B, file='figs/M4_B.png', width=16, height=16, units='cm')


# look at consistency within individuals `in_both_ppmi`

dat.bothppmi <- dat.all[ID %in% in_both_ppmi]

dat.bothppmi.wide <- dcast(dat.bothppmi, ID+UniProt+mode+sampletype+month+ARM+M4~STUDY, value.var='value')

ggplot(dat.bothppmi.wide[sampletype=='CSF' & month==0], aes(x=`PPMI1`, y=`PPMI2`)) +
geom_point() +
facet_wrap(~ID) +
theme_few(12) +
scale_x_continuous(breaks=0, labels=NULL) +
scale_y_continuous(breaks=0, labels=NULL) +
labs(title='CSF measurement reproducibility at baseline')

# Baseline correlation
correlations <- dat.bothppmi.wide[, list('cor'=cor(PPMI1, PPMI2, use='pairwise.complete.obs')), by=list(ID,sampletype,ARM,month)]
correlations[, 'cor_label' := substr(cor, 1, 4)]
dat.bothppmi.wide[, ID := factor(ID)]

plot_correlation <- function(i.sampletype, i.month) {
    g <- ggplot(dat.bothppmi.wide[sampletype==i.sampletype & month==i.month], aes(x=`PPMI1`, y=`PPMI2`, color=ARM)) +
    geom_point() +
    geom_text_repel(data=correlations[sampletype==i.sampletype & month==i.month], 
                color='black', 
                aes(x=Inf, hjust=1,
                    y=-Inf, vjust=0,
                    label=cor_label
                )
            ) +
    facet_wrap(~ID, drop=FALSE) +
    theme_few(10) +
    scale_x_continuous(breaks=0, labels=NULL) +
    scale_y_continuous(breaks=0, labels=NULL) +
    labs(title=paste0(i.sampletype, ' measurement reproducibility at month ', i.month)) +
    geom_abline(slope=1, intercept=0, linetype='dashed', alpha=0.3)
    ggsave(g, file=paste0('figs/correlation_', i.sampletype, '_month', i.month, '.png'),
                width=20, height=20, units='cm')
}
plot_correlation('PLA', 0)
plot_correlation('PLA', 12)
plot_correlation('PLA', 24)
plot_correlation('CSF', 0)
plot_correlation('CSF', 12)
plot_correlation('CSF', 24)

g.all_correlation <- ggplot(correlations, aes(x=1, y=cor, color=ARM)) +
    geom_beeswarm(shape=21, alpha=0.8) +
    facet_grid(sampletype~month, switch='x') +
    theme_few() +
    scale_x_continuous(breaks=0, labels=NULL) +
    ylim(0,1) +
    labs(x=NULL, y='Correlation')

ggsave(g.all_correlation, file='figs/correlation_all.png', width=25, height=15, units='cm')



## look at replication
# add replicate column per combination which should be 'unique'
# dat[, replicate := 1:.N, by=list(UniProt, participant_id, month, sampletype)]
# dat[UniProt=='P01375' & participant_id=='3004' & sampletype=='CSF' & month==0]
#    UniProt   value            mode PD_or_PP participant_id sampletype D month
# 1:  P01375 -1.2780 cardiometabolic       PP           3004        CSF 1     0
# 2:  P01375 -1.3305    inflammation       PP           3004        CSF 1     0
# 3:  P01375 -1.3890       neurology       PP           3004        CSF 1     0
# 4:  P01375 -0.7351        oncology       PP           3004        CSF 1     0
# 5:  P01375 -0.7945 cardiometabolic       PP           3004        CSF 2     0
# 6:  P01375 -0.8917    inflammation       PP           3004        CSF 2     0
# 7:  P01375 -0.5618       neurology       PP           3004        CSF 2     0
# 8:  P01375 -0.7578        oncology       PP           3004        CSF 2     0
# dat[, replicate := 1:.N, by=list(UniProt, participant_id, month, sampletype)]

## Super simple, take mean for every unique person/timepoint
dat.ag <- dat[, list('value'=mean(value, na.rm=T)), by=list(UniProt, participant_id, month, sampletype)]

# reduces rows from 4214380 to 3961229 (less 253150)

# transform to wide
dat.wide <- dcast(dat.ag, UniProt+participant_id+month~sampletype, value.var='value')


# load gene UNIPROT:SYMBOL table
gene_map_fn <- paste0('/data/CARD/projects/proteomicsXprogression_PD/data-allreleases/',
                        'OLINK-UNIPROT-GENE-MAPPING.txt')
gene_map <- fread(gene_map_fn, header=TRUE)

# ARM
# "Healthy Control"             "PD"
# ""                            "LBD"
# "Disease Control"             "Unknown"
# "Prodromal"                   "SWEDD"en
# "Genetic Registry PD"         "Genetic Cohort PD"
# "Genetic Registry Unaffected" "Genetic Cohort Unaffected"

# mis-coded?
# enrollment[PHENO==2 & ENROLL_STUDY_ARM=='Healthy Control']
#       ID BASELINE_AGE SEX PHENO ENROLL_STUDY_ARM COHORT
#  PP-3160           79   1     2  Healthy Control   PPMI
#  PP-3191           66   2     2  Healthy Control   PPMI
#  PP-3310           64   1     2  Healthy Control   PPMI
#  PP-3478           77   1     2  Healthy Control   PPMI

# LBD = lewy body dementia




enrollment[, c('group','ID') := tstrsplit(ID, split='-')]

# subset only participants in PPMI
ppmi_inds <- enrollment[COHORT == 'PPMI' & group=='PP']
setkey(ppmi_inds, ID)
# 425 PD
# 218 Healthy Control
setnames(dat, 'participant_id', 'ID')
setkey(dat, ID)



pdbp_inds <- enrollment[COHORT == 'PDBP' & group == 'PD']
setkey(pdbp_inds, ID)


ppmi <- merge(dat, ppmi_inds)
pdbp <- merge(dat, pdbp_inds)


# for PROTEIN DATA
# participant_id is either fully numeric (4 or 5 digits) or PD########"

#  enrollment[, .N, by=list(PHENO,ENROLL_STUDY_ARM)][order(PHENO, N)]
#     PHENO            ENROLL_STUDY_ARM    N
#  1:     0           Genetic Cohort PD    4
#  2:     0                       SWEDD    9
#  3:     0             Healthy Control   15
#  4:     0 Genetic Registry Unaffected   16
#  5:     0                          PD   27
#  6:     0   Genetic Cohort Unaffected   44
#  7:     0                   Prodromal   45
#  8:     0                               97
#  9:     0             Disease Control  155
# 10:     0                         LBD 2521

# 11:     1             Disease Control    1
# 12:     1           Genetic Cohort PD    1
# 13:     1                       SWEDD    2
# 14:     1                   Prodromal    4
# 15:     1 Genetic Registry Unaffected  226
# 16:     1   Genetic Cohort Unaffected  338
# 17:     1                              382
# 18:     1             Healthy Control 3358
# 19:     2                     Unknown    1
# 20:     2             Disease Control    3

# 21:     2 Genetic Registry Unaffected    3
# 22:     2             Healthy Control    4
# 23:     2                   Prodromal   15
# 24:     2   Genetic Cohort Unaffected   21
# 25:     2                       SWEDD   66
# 26:     2         Genetic Registry PD  196
# 27:     2           Genetic Cohort PD  267
# 28:     2                              273
# 29:     2                          PD 2678

enrollment[PHENO==2 & ENROLL_STUDY_ARM == 'PD']
enrollment[PHENO==1 & ENROLL_STUDY_ARM == 'Healthy Control']

# include rows with paired data
dat.wide <- dat.wide[!is.na(CSF) & ! is.na(PLA)]

library(ggplot2)
ggplot(dat.wide, aes(x=CSF, y=PLA)) +
    geom_point(alpha=0.2, shape=21) +
    facet_wrap(.~month)

