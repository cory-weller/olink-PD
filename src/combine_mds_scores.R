#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)

# data directory
dd='/data/CARD/projects/proteomicsXprogression_PD/data-allreleases/'

# merge mds1, 2, and 3
mds1 <- fread(paste0(dd, '/amp-pd-v25-clinical/clinical/MDS_UPDRS_Part_I.csv'))
mds1 <- mds1[, c('participant_id','visit_month','mds_updrs_part_i_summary_score')]
setkey(mds1, participant_id, visit_month)

mds2 <- fread(paste0(dd, '/amp-pd-v25-clinical/clinical/MDS_UPDRS_Part_II.csv'))
mds2 <- mds2[, c('participant_id','visit_month','mds_updrs_part_ii_summary_score')]
setkey(mds2, participant_id, visit_month)

mds3 <- fread(paste0(dd, '/amp-pd-v25-clinical/clinical/MDS_UPDRS_Part_III.csv'))
mds3 <- mds3[, c('participant_id','visit_month','mds_updrs_part_iii_summary_score')]
setkey(mds3, participant_id, visit_month)

merge1 <- merge(mds1, mds2, all=TRUE)
merge2 <- merge(merge1, mds3, all=TRUE)

setnames(merge2, 'mds_updrs_part_i_summary_score', 'mds_i')
setnames(merge2, 'mds_updrs_part_ii_summary_score', 'mds_ii')
setnames(merge2, 'mds_updrs_part_iii_summary_score', 'mds_iii')

dat <- melt(merge2, measure.vars=c('mds_i', 'mds_ii', 'mds_iii'),
            variable_name='mds_section', value.name='summary_score')

# some patients have multiple scores for a given time point, take average per month
dat <- dat[, list('summary_score'=mean(summary_score, na.rm=T)), by=list(participant_id, visit_month, variable)]
dat <- dat[! is.na(summary_score)]


# exclude visit months 0.5, -1 and -2
dat <- dat[! visit_month %in% c(-2, -1, 0.5)]

# output
fwrite(dat, file='reports/mds_scores.tsv', quote=F, row.names=F, col.names=T, sep='\t')