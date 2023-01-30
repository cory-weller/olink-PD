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

# calculate composite for those with info at all three months
hq_participants <- dat[! is.na(summary_score), .N, by=list(participant_id, visit_month)][N==3, participant_id]
dat.hq <- dat[participant_id %in% hq_participants]


dat.hq.total <- dat.hq[, list('total'=sum(summary_score)), by=list(participant_id, visit_month)]
dat.hq.total[, jittered := total + runif(.N, -0.5, 0.5)]

lowstart_pids <- unique(dat.hq.total[visit_month <= 6 & total <= 25, participant_id])

ggplot(dat.hq.total[participant_id %in% lowstart_pids],
    aes(x=visit_month, y=jittered, group=participant_id)) +
    geom_line(alpha=0.1) +
    labs(x='month', y='MDS composite score', title='score over time for initial score of 0')


###

# The following files relate to olink protein expression measured for participants

# './proteomics-CSF-PPEA-D02/olink-explore/protein-expression/CSF-PPEA-D02_matrix_cardiometabolic.csv'
# './proteomics-CSF-PPEA-D02/olink-explore/protein-expression/CSF-PPEA-D02_matrix_neurology.csv'
# './proteomics-CSF-PPEA-D02/olink-explore/protein-expression/CSF-PPEA-D02_matrix_inflammation.csv'
# './proteomics-CSF-PPEA-D02/olink-explore/protein-expression/CSF-PPEA-D02_matrix_oncology.csv'

# './proteomics-CSF-PPEA-D01/olink-explore/protein-expression/CSF-PPEA-D01_matrix_cardiometabolic.csv'
# './proteomics-CSF-PPEA-D01/olink-explore/protein-expression/CSF-PPEA-D01_matrix_inflammation.csv'
# './proteomics-CSF-PPEA-D01/olink-explore/protein-expression/CSF-PPEA-D01_matrix_oncology.csv'
# './proteomics-CSF-PPEA-D01/olink-explore/protein-expression/CSF-PPEA-D01_matrix_neurology.csv'

# './proteomics-PLA-PPEA-D02/olink-explore/protein-expression/PLA-PPEA-D02_matrix_neurology.csv'
# './proteomics-PLA-PPEA-D02/olink-explore/protein-expression/PLA-PPEA-D02_matrix_cardiometabolic.csv'
# './proteomics-PLA-PPEA-D02/olink-explore/protein-expression/PLA-PPEA-D02_matrix_inflammation.csv'
# './proteomics-PLA-PPEA-D02/olink-explore/protein-expression/PLA-PPEA-D02_matrix_oncology.csv'

# './proteomics-PLA-PPEA-D01/olink-explore/protein-expression/PLA-PPEA-D01_matrix_neurology.csv'
# './proteomics-PLA-PPEA-D01/olink-explore/protein-expression/PLA-PPEA-D01_matrix_cardiometabolic.csv'
# './proteomics-PLA-PPEA-D01/olink-explore/protein-expression/PLA-PPEA-D01_matrix_inflammation.csv'
